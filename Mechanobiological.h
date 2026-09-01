#pragma once

#include <FEBioMix/FEMultiphasic.h>
#include <FEBioMix/FESolutesMaterialPoint.h>
#include <FEBioMix/FESolute.h>
#include <FEBioMix/FEChemicalReaction.h>
#include <FEBioMix/FESoluteDiffusivity.h>

#include <FECore/FEModel.h>
#include <FECore/FEMaterialPoint.h>
#include <FECore/FELoadController.h>
#include <FECore/FEModelParam.h>

#include <vector>


/**
 * @brief Material-point data used by the mechanobiological material.
 *
 * Stores history variables required to update the referential fixed charge
 * density during long-term tissue remodeling.
 */

class FEMechanobiological;


class FEBIOMIX_API FEMechanobiologicalPoint : public FEMaterialPointData
{
public:
	FEMechanobiologicalPoint(FEMechanobiological* pm, FEMaterialPointData* pt)
		: FEMaterialPointData(pt), m_pMat(pm) {}

	FEMaterialPointData* Copy();

	void Serialize(DumpStream& ar);

	void Init();


public:
	FEMechanobiological* m_pMat = nullptr;
public:	

	// History variables associated with the mechanobiological evolution.
	vector<mat3d> m_Fi;
	vector<double> m_Ji;

	// Referential fixed charge density.
	double FCDr = 0.0;    

	// Cumulative change in referential fixed charge density.
	double FCD_loss = 0.0; 

	// Time step at which the FCD was last updated.
	int lastTimeStep = -1;  

};


// ============================================================================
// Mechanobiological multiphasic material
// ============================================================================


/**
 * @brief Multiphasic material with long-term mechanobiological remodeling.
 *
 * Extends FEBio's multiphasic material by coupling the evolution of the
 * referential fixed charge density to extracellular matrix biosynthesis
 * and degradation.
 *
 * The referential FCD evolves according to
 *
 *     d(FCD_r)/dt = -S_bio + D_deg
 *
 * where S_bio is the biosynthesis contribution and D_deg is the degradation
 * contribution.
 * 
 * @note The degeneration grade is represented as a global IVD-level
 * variable rather than as an independently evolving local field. The global
 * grade is determined from the most advanced local degeneration state in the
 * nucleus pulposus (NP), based on the local referential fixed charge density
 * (or equivalently the local proteoglycan content).
 *
 * In the numerical implementation, the grade is introduced as an auxiliary
 * solute with effectively infinite diffusivity. Its value is prescribed
 * through a load controller, which evaluates the local NP degeneration state
 * and imposes the resulting global grade throughout the IVD. This allows
 * the same degeneration-dependent constitutive laws to be used consistently
 * across all tissue regions without solving an additional transport problem
 * for the grade.
 */


class FEMechanobiological : public FEMultiphasic
{
public:

	FEMechanobiological(FEModel* pfem);

	FEMaterialPointData* CreateMaterialPointData() override;

	void UpdateSolidBoundMolecules(FEMaterialPoint& mp) override;

	double FixedChargeDensity(FEMaterialPoint& mp);

	tens4ds Tangent(FEMaterialPoint& mp);

	void PartitionCoefficientFunctions(FEMaterialPoint& mp, vector<double>& kappa,
		vector<double>& dkdJ,
		vector< vector<double> >& dkdc,
		vector< vector<double> >& dkdr,
		vector< vector<double> >& dkdJr,
		vector< vector< vector<double> > >& dkdrc);

	double Degradation(FEMaterialPoint& mp);
	double Biosynthesis(FEMaterialPoint& mp);
	double GetReferentialFixedChargeDensity(const FEMaterialPoint& mp);

public:
	int m_region; //!< Tissue region: 0 = nucleus pulposus, 1 = annulus fibrosus, 2 = cartilage endplate.
	double m_FCD0;   //!< Initial referential fixed charge density.
	double m_delta;  //!< Degeneration-dependent reduction factor for biosynthesis.
	DECLARE_FECORE_CLASS();

};

// ============================================================================
// Degeneration grade controller
// ============================================================================

/**
 * @brief Load controller returning the global degeneration grade.
 *
 * The degeneration grade is determined from the maximum local reduction
 * in referential fixed charge density (FCD) within the nucleus pulposus.
 *
 * The controller is used to convert the local mechanobiological state of
 * the nucleus pulposus into a global degeneration grade that can be used
 * as a model output or as a loading/control variable.
 */
class Grade : public FELoadController
{
public:
	//! Constructor
	Grade(FEModel* fem);

	//! Returns the current global degeneration grade.
	double GetValue(double time) override;

protected:
	//! Initial referential fixed charge density used as the reference state.
	double m_FCD0;

	DECLARE_FECORE_CLASS();
};


// ============================================================================
// Effective OCF stiffness controller
// ============================================================================

/**
 * @brief Load controller returning the effective stiffness of oriented
 *        collagen fibers (OCF).
 *
 * The effective OCF stiffness is computed as a function of the current
 * degeneration state and the corresponding stiffness parameter.
 */
class OCF_stiffness_eff : public FELoadController
{
public:
	//! Constructor
	OCF_stiffness_eff(FEModel* fem);

	//! Returns the current effective OCF stiffness.
	double GetValue(double time) override;

protected:
	//! Reference lamellar stiffness parameter.
	double m_lam;

	DECLARE_FECORE_CLASS();
};

/**
 * @brief Aerobic glucose/oxygen metabolism reaction.
 *
 * Computes the oxygen consumption associated with aerobic cellular metabolism
 * as a function of local oxygen and glucose concentrations, cell density,
 * tissue region, and degeneration grade.
 *
 * The analytical derivatives required by FEBio's Newton solver are provided
 * with respect to deformation and solute concentrations.
 */


 // ============================================================================
 // Aerobic metabolism reaction
 // ============================================================================
/**
 * @brief Aerobic oxygen consumption reaction.
 *
 * Computes oxygen consumption as a function of:
 * - local glucose concentration,
 * - local oxygen concentration,
 * - effective cell density,
 * - tissue region,
 * - degeneration grade.
 *
 * Analytical tangents are provided with respect to deformation and effective
 * solute concentrations.
 */
class FEBIOMIX_API Aerobic_Reaction : public FEChemicalReaction
{
public:
	Aerobic_Reaction(FEModel* pfem);
	bool Init() override;
	double ReactionSupply(FEMaterialPoint& pt) override;
	mat3ds Tangent_ReactionSupply_Strain(FEMaterialPoint& pt) override;
	double Tangent_ReactionSupply_Pressure(FEMaterialPoint& pt) override;
	double Tangent_ReactionSupply_Concentration(FEMaterialPoint& pt, const int sol) override;

public:
	int		m_Rid;			//!< local id of reactant
	int		m_Pid;			//!< local id of product
	bool	m_Rtype;		//!< flag for reactant type (solute = false, sbm = true)
	int m_idGlucose;
	int m_idOxygen;
	int m_region;// int: 0=NP, 1=AF, 2=CEP;

	DECLARE_FECORE_CLASS();
};


// ============================================================================
// Anaerobic metabolism reaction
// ============================================================================

/**
 * @brief Anaerobic glucose consumption reaction.
 *
 * Computes glucose consumption associated with anaerobic metabolism.
 *
 * The model includes:
 * - basal lactate production regulated by glucose and oxygen availability,
 * - oxygen-dependent metabolic partitioning,
 * - degeneration-dependent maintenance metabolism,
 * - glucose-dependent effective cell density.
 *
 * The glucose consumption rate is obtained from the lactate production rate
 * using the stoichiometric relation
 *
 *     Q_glucose = Q_lactate / 2.
 *
 * Analytical tangents are provided with respect to deformation and effective
 * solute concentrations.
 */

class FEBIOMIX_API Anaerobic_Reaction : public FEChemicalReaction
{
public:
	Anaerobic_Reaction(FEModel* pfem);
	bool Init() override;
	double ReactionSupply(FEMaterialPoint& pt) override;
	mat3ds Tangent_ReactionSupply_Strain(FEMaterialPoint& pt) override;
	double Tangent_ReactionSupply_Pressure(FEMaterialPoint& pt) override;
	double Tangent_ReactionSupply_Concentration(FEMaterialPoint& pt, const int sol) override;

public:
	int		m_Rid;			//!< local id of reactant
	int		m_Pid;			//!< local id of product
	bool	m_Rtype;		//!< flag for reactant type (solute = false, sbm = true)
	int		m_region;		// int: 0=NP, 1=AF, 2=CEP
	int m_idGlucose;
	int m_idOxygen;

	DECLARE_FECORE_CLASS();
};

// ============================================================================
// Degeneration-dependent solute diffusivity
// ============================================================================

/**
 * @brief Isotropic solute diffusivity modified by degeneration.
 *
 * The diffusivity is reduced as a function of the local degeneration grade,
 * which is used here as a surrogate for cartilage endplate calcification.
 */

class FEBIOMIX_API FEDiffusivity_Calcification : public FESoluteDiffusivity
{
public:
	FEDiffusivity_Calcification(FEModel* pfem);

	double Free_Diffusivity(FEMaterialPoint& pt) override;

	double Tangent_Free_Diffusivity_Concentration(FEMaterialPoint& mp, const int isol) override;

	mat3ds Diffusivity(FEMaterialPoint& pt) override;

	tens4dmm Tangent_Diffusivity_Strain(FEMaterialPoint& mp) override;

	mat3ds Tangent_Diffusivity_Concentration(FEMaterialPoint& mp, const int isol) override;

	bool Validate() override;

public:
	FEParamDouble	m_free_diff;	//!< free diffusivity
	FEParamDouble	m_diff;			//!< diffusivity

	DECLARE_FECORE_CLASS();
};


