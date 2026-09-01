#include "Mechanobiological.h"

#include "FEMechanobiological.h"

#include <FEBioMix/FESolutesMaterialPoint.h>
#include <FEBioMix/FESoluteDiffusivity.h>

#include <FECore/FEElement.h>
#include <FECore/FEDomain.h>
#include <FECore/FEMesh.h>
#include <FECore/FEMaterial.h>
#include <FECore/FEMaterialPoint.h>
#include <FECore/FEShellElement.h>
#include <FECore/FELoadController.h>
#include <FECore/log.h>

#include <algorithm>
#include <cmath>

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

// ============================================================================
// Mechanobiological multiphasic material
// ============================================================================

BEGIN_FECORE_CLASS(FEMechanobiological, FEMultiphasic)
    ADD_PARAMETER(m_region, "region");
    ADD_PARAMETER(m_FCD0, "cFr_init");
    ADD_PARAMETER(m_delta, "delta");
END_FECORE_CLASS();


FEMechanobiological::FEMechanobiological(FEModel* pfem) : FEMultiphasic(pfem)
{
}

// ============================================================================
// Mechanobiological material point
// ============================================================================

FEMaterialPointData* FEMechanobiological::CreateMaterialPointData()
{
    return new FEMechanobiologicalPoint
           (this, new FESolutesMaterialPoint
           (new FEBiphasicMaterialPoint
           (m_pSolid->CreateMaterialPointData())
            )
            );
}


void FEMechanobiological::UpdateSolidBoundMolecules(FEMaterialPoint& mp)
{
    const FETimeInfo& ti = GetFEModel()->GetTime();
    double dt = CurrentTimeIncrement();
    double t = ti.currentTime;

    int nreact = (int)Reactions();
    int mreact = (int)MembraneReactions();
    if (nreact || mreact) {
        // for chemical reactions involving solid-bound molecules,
        // update their concentration
        // multiphasic material point data
        FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
        FEBiphasicMaterialPoint& ppt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
        FESolutesMaterialPoint& spt = *(mp.ExtractData<FESolutesMaterialPoint>());
        FEShellElement* sel = dynamic_cast<FEShellElement*>(mp.m_elem);
        double h = (sel) ? sel->Evaluate(sel->m_ht, mp.m_index) : 0;  

        double phi0 = ppt.m_phi0t;
        int nsbm = SBMs();
        int nsol = Solutes();
        // create a temporary container for spt.m_sbmr so that this variable remains
        // unchanged as the reaction rates get calculated for each sbm.
        vector<double> sbmr = spt.m_sbmr;
        for (int isbm = 0; isbm < nsbm; ++isbm) {
            spt.m_sbmrhat[isbm] = 0;
            // combine the molar supplies from all the reactions
            for (int k = 0; k < nreact; ++k) {
                double zetahat = GetReaction(k)->ReactionSupply(mp);
                double v = GetReaction(k)->m_v[nsol + isbm];
                // remember to convert from molar supply to referential mass supply
                spt.m_sbmrhat[isbm] += (pt.m_J - phi0) * SBMMolarMass(isbm) * v * zetahat;
            }
            for (int k = 0; k < mreact; ++k) {
                double zetahat = GetMembraneReaction(k)->ReactionSupply(mp);
                double v = GetMembraneReaction(k)->m_v[nsol + isbm];
                // remember to convert from molar supply to referential mass supply
                spt.m_sbmrhat[isbm] += pt.m_J / h * SBMMolarMass(isbm) * v * zetahat;
            }
            // perform the time integration (midpoint rule)
            sbmr[isbm] = spt.m_sbmrp[isbm] + dt * (spt.m_sbmrhat[isbm] + spt.m_sbmrhatp[isbm]) / 2;
            // check bounds
            if (sbmr[isbm] < spt.m_sbmrmin[isbm])
                sbmr[isbm] = spt.m_sbmrmin[isbm];
            if ((spt.m_sbmrmax[isbm] > 0) && (sbmr[isbm] > spt.m_sbmrmax[isbm]))
                sbmr[isbm] = spt.m_sbmrmax[isbm];
        }
        // now update spt.m_sbmr
        spt.m_sbmr = sbmr;
    }
}

FEMaterialPointData* FEMechanobiologicalPoint::Copy()
{
    FEMechanobiologicalPoint* pt = new FEMechanobiologicalPoint(*this);
    if (m_pNext) pt->m_pNext = m_pNext->Copy();
    return pt;
}

void FEMechanobiologicalPoint::Serialize(DumpStream& ar)
{
    FEMaterialPointData::Serialize(ar);
}

void FEMechanobiologicalPoint::Init()
{
    FEMaterialPointData::Init();

    m_Fi.clear();
    m_Ji.clear();

    FESolutesMaterialPoint& spt = *((*this).ExtractData<FESolutesMaterialPoint>());
    mat3d Fi(mat3dd(1.0));
    double Ji = 1;
    m_Fi.push_back(Fi);
    m_Ji.push_back(Ji);
}

// ============================================================================
// Fixed charge density
// ============================================================================


double FEMechanobiological::FixedChargeDensity(FEMaterialPoint& pt)
{
    FEElasticMaterialPoint& et = *pt.ExtractData<FEElasticMaterialPoint>();
    FEBiphasicMaterialPoint& bt = *pt.ExtractData<FEBiphasicMaterialPoint>();
    FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();

    auto& mpd = *pt.ExtractData<FEMechanobiologicalPoint>();

    const FETimeInfo& ti = GetFEModel()->GetTime();
    double dt = ti.timeIncrement ; 
    double J = et.m_J;
    double phi0 = bt.m_phi0t;
    double ce = 0;

    // add contribution from charged solid-bound molecules
    for (int isbm = 0; isbm < (int)m_pSBM.size(); ++isbm)
        ce += SBMChargeNumber(isbm) * spt.m_sbmr[isbm] / SBMMolarMass(isbm);

    // Update the referential FCD once per time step.
    double Delta_FCD = 0;
    if (ti.timeStep != mpd.lastTimeStep)
    {
        Delta_FCD = (-Biosynthesis(pt) + Degradation(pt)) * dt / (3600.0 * 24.0 * 365.0); //s -> years
        mpd.FCD_loss += Delta_FCD;
        mpd.FCDr = m_cFr(pt) + mpd.FCD_loss;
        mpd.lastTimeStep = ti.timeStep;
    }
    double cF = (mpd.FCDr * (1 - bt.m_phi0) + ce) / (J - bt.m_phi0); 

    // add contribution from solid-bound 'solutes'
    const int nsol = (int)m_pSolute.size();
    for (int isol = 0; isol < nsol; ++isol)
        if (spt.m_bsb[isol]) cF += spt.m_ca[isol] * m_pSolute[isol]->ChargeNumber();

    return cF;
}

// ============================================================================
// Degradation
// ============================================================================

double FEMechanobiological::Degradation(FEMaterialPoint& pt)
{
    const FETimeInfo& ti = GetFEModel()->GetTime();
    double t = ti.currentTime / (365 * 24 * 3600); // s-->years, t=A-A0

    auto& mpd = *pt.ExtractData<FEMechanobiologicalPoint>();

    FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
    double G0 = 1.5; //initial grade (representing healthy and young state)
    double G = std::min(5.0, std::max(G0, (spt.m_c.empty() ? G0 : spt.m_c[3]))); //current grade

    double kA0 = 0.17;
    double kAinf = 0.12;
    double kPG_G = 0.06;
    double tau_A = 20.0;
    double tau_G = 40.0;

    double fA = kAinf / kA0 + (1 - kAinf / kA0) * exp(-t/ tau_A); //Equation A.7
    double kPG_ND = kA0 * fA; //Equation A.6

    double fG = (G - G0) * exp(-t / tau_G); //Equation A.9
    double deltakPG = kPG_G * fG; //Equation A.8

    double kPG = kPG_ND + deltakPG; //
    double kFCD = kPG;

    double D_FCD = -kFCD * mpd.FCDr; //Equation 25

    return D_FCD;
}

// ============================================================================
// Biosynthesis
// ============================================================================

double FEMechanobiological::Biosynthesis(FEMaterialPoint& mp)
{
    FESolutesMaterialPoint& spt = *(mp.ExtractData<FESolutesMaterialPoint>());
    double cG = (spt.m_c.empty() ? 0.0 : spt.m_c[0]); //glucose concentration
    double cO = (spt.m_c.empty() ? 0.0 : spt.m_c[2]); //oxygen concentration
    double G0 = 1.5; //initial grade
    double G = std::min(5.0, std::max(G0, (spt.m_c.empty() ? G0 : spt.m_c[3])));
     
    const FETimeInfo& ti = GetFEModel()->GetTime();
    double t = ti.currentTime / (3600 * 24 * 365); //s-->years

    double kA0 = 0.17;
    double kAinf = 0.12;
    double tau_A = 20.0;

    double fA = kAinf / kA0 + (1 - kAinf / kA0) * exp(-t / tau_A); //Equation A.7

    double Ylac = 1.0;
    double YO2 = 5.0;
    double cG0 = 25.0;
    double cg_crit = 1.0; //glucose threshold
    double QO2_max = 18.0;
    double KO2 = 0.0127;
    double Qlac_max = 0.0;
    double Ap = 0.0;
    double Bp = 0.0;
    double alpha_p = 0.0;
    double Beta_p = 0.0;

    double rho_cell_eff = 0.0;
    double rho_cell0 = 0.0; 

    if (m_region == 0) { //NP
        rho_cell_eff = std::min(2000.0, 2000.0 * cG / cg_crit);
        rho_cell0 = 2000.0;
        Qlac_max = 160;
        Ap = 1.31;
        Bp = 0.12;
        alpha_p = 35;
        Beta_p = 0.07;
    }
    else if (m_region == 1) { //AF
        rho_cell_eff = std::min(2000.0, 2000.0 * cG / cg_crit);
        rho_cell0 = 2000.0;
        Qlac_max = 40;
        Ap = 2.98;
        Bp = 0.28;
        alpha_p = 8;
        Beta_p = 0.02;
    }
    else if (m_region == 2) { //CEP
        rho_cell_eff = std::min(10000.0, 10000.0 * cG / cg_crit);
        rho_cell0 = 10000.0;
        Qlac_max = 40;
        Ap = 2.98;
        Bp = 0.28;
        alpha_p = 8;
        Beta_p = 0.02;
    }

    double QO2_bio = QO2_max * cO / (KO2 + cO);
    double delta_Qlac = 5 * (QO2_max - QO2_bio);
    double f_nut = 1 - Ap * exp(-cG * cG / alpha_p) * exp(-cO / Beta_p) - Bp * (cG0 - cG) / cG0;
    double Qlac_bio = Qlac_max * f_nut + delta_Qlac;
    double factor = pow(10, -6) / 3600.0; // Conversion from micromolar/h to mol/m^3/s, if applicable to the parameterization used in the model.
    double RATP_bio = (YO2 * QO2_bio + Ylac * Qlac_bio)*factor; //Equation 24
    double RATP_bio0= (YO2 * QO2_max + Ylac * Qlac_max) * factor;
    double fG = 1 / (1 + m_delta * (G - G0));
    double eta0 = -m_FCD0 * kA0 / (RATP_bio0 * rho_cell0);
    double eta_PG = eta0 * fA * fG;

    double S_FCD = eta_PG * rho_cell_eff * RATP_bio; //Equation 26

    return S_FCD;
}

double FEMechanobiological::GetReferentialFixedChargeDensity(const FEMaterialPoint& pt)
{
    auto& mpd = *pt.ExtractData<FEMechanobiologicalPoint>();

    return mpd.FCDr;
}

// ============================================================================
// Degeneration grade controller
// ============================================================================

BEGIN_FECORE_CLASS(Grade, FELoadController);
ADD_PARAMETER(m_FCD0, "m_FCDR0");
END_FECORE_CLASS();


Grade::Grade(FEModel* fem)
    : FELoadController(fem)
{
}
 
//-----------------------------------------------------------------------------
double Grade::GetValue(double time)
{
    FEMesh& mesh = GetFEModel()->GetMesh();
    double G0 = 1.5;
    double G = G0;

    // Search for the maximum local degeneration grade in the NP.
    for (int d = 0; d < mesh.Domains(); ++d)
    {
        FEDomain& dom = mesh.Domain(d);
        FEMaterial* mat = dom.GetMaterial();

        auto* myMat = dynamic_cast<FEMechanobiological*>(mat);
        if (!myMat) continue;
        if (myMat->m_region != 0) continue; // NP only

        for (int i = 0; i < dom.Elements(); ++i)
        {
            FEElement& el = dom.ElementRef(i);
            int nGauss = el.GaussPoints();

            for (int j = 0; j < nGauss; ++j)
            {
                FEMaterialPoint* mp = el.GetMaterialPoint(j);
                if (!mp) continue;

                auto* mpd = mp->ExtractData<FEMechanobiologicalPoint>();
                if (!mpd) continue;
                double fcdr = mpd->FCDr;
                double ratio = fcdr / m_FCD0;
                if (ratio > 0.0)
                {
                    double G_local = -1.692 * log(ratio) + G0;
                    G_local = std::min(5.0, std::max(G0, G_local));
                    G = std::max(G, G_local);
                }

            }
        }
    }

    
    return G;
}

// ============================================================================
// OCF effective stiffness controller
// ============================================================================


BEGIN_FECORE_CLASS(OCF_stiffness_eff, FELoadController)
ADD_PARAMETER(m_lam, "lam");
END_FECORE_CLASS();

OCF_stiffness_eff::OCF_stiffness_eff(FEModel* fem)
    : FELoadController(fem)
{
}

//-----------------------------------------------------------------------------
double OCF_stiffness_eff::GetValue(double time)
{
    FEMesh& mesh = GetFEModel()->GetMesh();
    double G = 0.0;

    int d = 0; 
    FEDomain& dom = mesh.Domain(d);
    FEMaterial* mat = dom.GetMaterial();
    auto* myMat = dynamic_cast<FEMechanobiological*>(mat);

    int i = 0;
    FEElement& el = dom.ElementRef(i);
    int nGauss = el.GaussPoints();

    int j = 0; 
    FEMaterialPoint* mp = el.GetMaterialPoint(j);
    FESolutesMaterialPoint* spts = mp->ExtractData<FESolutesMaterialPoint>();
    G = spts->m_ca[3]; //grade
    double G1 = 6.0;
    double G0 = 1.5;
    double g_theta = 1.0 - std::exp(G - G0) / std::exp(G1 - G0); //Equation A.20, decrease in stiffness with degeneration
    double Einner = 20.98; //!< Healthy effective OCF stiffness in the inner AF.
    double Eouter = 79.57; //!< Healthy effective OCF stiffness in the outer AF.
    double n_lam = 15.0;   //!< Total number of AF lamellae.
    double E_OCF = (Einner +  (Eouter - Einner)* (m_lam - 1) / (n_lam -1)) * g_theta;

    if (!std::isfinite(E_OCF) || E_OCF <= 0.0) return 1.0;

    return E_OCF;
}

// ============================================================================
// Aerobic reaction
// ============================================================================

BEGIN_FECORE_CLASS(Aerobic_Reaction, FEChemicalReaction)
ADD_PARAMETER(m_idGlucose, "id_glucose");
ADD_PARAMETER(m_idOxygen, "id_oxygen");
ADD_PARAMETER(m_region, "IVD region");
END_FECORE_CLASS();

Aerobic_Reaction::Aerobic_Reaction(FEModel* pfem) : FEChemicalReaction(pfem)
{
    m_Rid = m_Pid = -1;
    m_Rtype = false;
}

bool Aerobic_Reaction::Init()
{
    if (FEChemicalReaction::Init() == false) return false;
    if (m_solR.size() + m_sbmR.size() > 1) {
        feLogError("Provide only one vR for this reaction");
        return false;
    }

    if (m_solP.size() + m_sbmP.size() > 1) {
        feLogError("Provide only one vP for this reaction");
        return false;
    }

    const int ntot = (int)m_v.size();
    for (int itot = 0; itot < ntot; itot++) {
        if (m_vR[itot] > 0) m_Rid = itot;
        if (m_vP[itot] > 0) m_Pid = itot;
    }

    if (m_Rid == -1) {
        feLogError("Provide vR for the reactant");
        return false;
    }

    if (m_Rid >= m_nsol) m_Rtype = true;

    return true;
}

double Aerobic_Reaction::ReactionSupply(FEMaterialPoint& pt)
{
    // Nutritional environment
    FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
    double cG = m_psm->GetActualSoluteConcentration(pt, m_idGlucose); //glucose concentration
    double cO = m_psm->GetActualSoluteConcentration(pt, m_idOxygen); //oxygen concentration

    // Degeneration grade
    double G0 = 1.5; //initial grade
    double G = std::min(5.0, std::max(G0, (spt.m_c.empty() ? G0 : spt.m_c[3]))); //current grade

    // Cell density
    double cg_crit = 1.0; //critical glucose threshold triggering cell death
    double rho_cell_eff = 0.0; //effective cell density
    if (m_region == 0) { //NP
        rho_cell_eff = std::min(2000.0, 2000.0 * cG / cg_crit);
    }
    else if (m_region == 1) { //AF
        rho_cell_eff = std::min(2000.0, 2000.0 * cG / cg_crit);
    }
    else if (m_region == 2) { //CEP
        rho_cell_eff = std::min(10000.0, 10000.0 * cG / cg_crit);
    }

    // Parameters
    double QO2_max = 18.0;
    double KO2 = 0.0127;
    double alpha_m = 2.0;
    double beta_m = 0.06;
    double delta_m = 0.6;
    double YO2 = 5.0;
    double Am = 152.5; 

    // Basal oxygen consumption (per cell)
    double QO2_bio = QO2_max * cO / (KO2 + cO);

    // Maintenance oxygen consumption (per cell)
    double QO2_maint = Am * delta_m / YO2 * (1.0 + exp(-cG / alpha_m)) * (1.0 - exp(-cO / beta_m)) * (G - G0); //=0 if G=G0

    // Total oxygen consumption
    double factor = pow(10, -6) / 3600.0; //unit conversion
    double QO2 = (QO2_bio + QO2_maint)* rho_cell_eff*factor;

    return QO2;
}

// ============================================================================
// Aerobic reaction tangent: deformation
// ============================================================================


mat3ds Aerobic_Reaction::Tangent_ReactionSupply_Strain(FEMaterialPoint& pt)
{
    double cG = m_psm->GetActualSoluteConcentration(pt, m_idGlucose);
    double cO = m_psm->GetActualSoluteConcentration(pt, m_idOxygen);
    double cGe = m_psm->GetEffectiveSoluteConcentration(pt, m_idGlucose);
    double cOe = m_psm->GetEffectiveSoluteConcentration(pt, m_idOxygen);

    double dcGdJ = m_psm->dkdJ(pt, m_idGlucose) * cGe;
    double dcOdJ = m_psm->dkdJ(pt, m_idOxygen) * cOe;

    FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
    double G0 = 1.5;
    double G = std::min(5.0, std::max(G0, (spt.m_c.empty() ? G0 : spt.m_c[3])));

    // Effective cell density and its derivative
    double rho_cell_eff = 0.0;
    double cg_crit = 1.0;
    double drho_dcG = 0.0;
    double rho_cell0 = 0.0;

    if (m_region == 0) { //NP
        rho_cell0 = 2000.0;
        rho_cell_eff = std::min(2000.0, 2000.0 * cG / cg_crit);
    }
    else if (m_region == 1) { //AF
        rho_cell0 = 2000.0;
        rho_cell_eff = std::min(2000.0, 2000.0 * cG / cg_crit);
    }
    else if (m_region == 2) { //CEP
        rho_cell0 = 10000.0;
        rho_cell_eff = std::min(10000.0, 10000.0 * cG / cg_crit);
    }

    if (cG < cg_crit)
    {
        drho_dcG = rho_cell0 / cg_crit;
    }

    // Parameters
    double QO2_max = 18.0;
    double KO2 = 0.0127;
    double alpha_m = 2.0;
    double beta_m = 0.06;
    double delta_m = 0.6;
    double YO2 = 5;
    double Am = 152.5;

    // Basal oxygen consumption
    double QO2_bio = QO2_max * cO / (KO2 + cO) * rho_cell_eff;
    double dQO2bio_dcO = QO2_max * KO2 / pow(KO2 + cO, 2) * rho_cell_eff;
    double dQO2bio_dcG = QO2_max * cO / (KO2 + cO) * drho_dcG;

    // Maintenance oxygen consumption
    double expG = exp(-cG / alpha_m);
    double expO = exp(-cO / beta_m);
    double C = Am * delta_m / YO2 * (G - G0);
    double QO2_maint = C * (1.0 + expG) * (1.0 - expO) * rho_cell_eff;

    double dQO2maint_dcG = -C * expG / alpha_m * (1.0 - expO) * rho_cell_eff + C * (1.0 + expG) * (1.0 - expO) * drho_dcG;
    double dQO2maint_dcO = C * (1.0 + expG) * expO / beta_m * rho_cell_eff;

    // Total derivatives with respect to actual concentrations
    double dQ_dcG = dQO2bio_dcG + dQO2maint_dcG;
    double dQ_dcO = dQO2bio_dcO + dQO2maint_dcO;

    // Chain rule with respect to J
    double factor = pow(10, -6) / 3600.0; //unit conversion
    double dQ_dJ = (dQ_dcG * dcGdJ + dQ_dcO * dcOdJ) * factor;

    return mat3dd(1.0) * dQ_dJ;
}

// ============================================================================
// Aerobic reaction tangent: pressure
// ============================================================================

double Aerobic_Reaction::Tangent_ReactionSupply_Pressure(FEMaterialPoint& pt)
{
    return 0;
}

// ============================================================================
// Aerobic reaction tangent: effective solute concentration
// ============================================================================


double Aerobic_Reaction::Tangent_ReactionSupply_Concentration(FEMaterialPoint& pt, const int sol)
{
    double cG = m_psm->GetActualSoluteConcentration(pt, m_idGlucose);
    double cO = m_psm->GetActualSoluteConcentration(pt, m_idOxygen);
    double cGe = m_psm->GetEffectiveSoluteConcentration(pt, m_idGlucose);
    double cOe = m_psm->GetEffectiveSoluteConcentration(pt, m_idOxygen);

    // Partition coefficients
    double kG = m_psm->GetPartitionCoefficient(pt, m_idGlucose);
    double kO = m_psm->GetPartitionCoefficient(pt, m_idOxygen);
    double dkGdc = m_psm->dkdc(pt, m_idGlucose, m_idGlucose);
    double dkOdc = m_psm->dkdc(pt, m_idOxygen, m_idOxygen);

    FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();

    double G0 = 1.5; 
    double G = std::min(5.0, std::max(G0, (spt.m_c.empty() ? G0 : spt.m_c[3])));

    // Effective cell density and its derivative
    double rho_cell_eff = 0.0;
    double cg_crit = 1.0;
    double drho_dcG = 0.0;
    double rho_cell0 = 0.0;

    if (m_region == 0) { //NP
        rho_cell0 = 2000.0;
        rho_cell_eff = std::min(2000.0, 2000.0 * cG / cg_crit);
    }
    else if (m_region == 1) { //AF
        rho_cell0 = 2000.0;
        rho_cell_eff = std::min(2000.0, 2000.0 * cG / cg_crit);
    }
    else if (m_region == 2) { //CEP
        rho_cell0 = 10000.0;
        rho_cell_eff = std::min(10000.0, 10000.0 * cG / cg_crit);
    }

    if (cG < cg_crit)
    {
        drho_dcG = rho_cell0 / cg_crit;
    }

    // Parameters
    double QO2_max = 18.0;
    double KO2 = 0.0127;
    double alpha_m = 2.0;
    double beta_m = 0.06;
    double delta_m = 0.6;
    double YO2 = 5.0;
    double Am = 152.5;

    // Basal oxygen consumption
    double QO2_bio = QO2_max * cO / (KO2 + cO) * rho_cell_eff;
    double dQO2bio_dcO = QO2_max * KO2 / pow(KO2 + cO, 2) * rho_cell_eff;
    double dQO2bio_dcG = QO2_max * cO / (KO2 + cO) * drho_dcG;

    // Maintenance oxygen consumption
    double expG = exp(-cG / alpha_m);
    double expO = exp(-cO / beta_m);
    double C = Am * delta_m / YO2 * (G - G0);
    double QO2_maint = C * (1.0 + expG) * (1.0 - expO) * rho_cell_eff;

    double dQO2maint_dcG = -C * expG / alpha_m * (1.0 - expO) * rho_cell_eff + C * (1.0 + expG) * (1.0 - expO) * drho_dcG;
    double dQO2maint_dcO = C * (1.0 + expG) * expO / beta_m * rho_cell_eff;

    // Total derivatives with respect to actual concentrations
    double dQ_dcG = dQO2bio_dcG + dQO2maint_dcG;
    double dQ_dcO = dQO2bio_dcO + dQO2maint_dcO;

    double factor = pow(10, -6) / 3600.0; //unit conversion

    // Derivative with respect to glucose
    if (sol == m_idGlucose)
    {
        return dQ_dcG * (kG + cGe * dkGdc) * factor;
    }

    // Derivative with respect to oxygen
    if (sol == m_idOxygen)
    {
        return dQ_dcO * (kO + cOe * dkOdc) * factor;
    }

    return 0.0;
}

// ============================================================================
// Anaerobic reaction
// ============================================================================

BEGIN_FECORE_CLASS(Anaerobic_Reaction, FEChemicalReaction)
ADD_PARAMETER(m_region, "region");// int: 0=NP, 1=AF ou CEP
ADD_PARAMETER(m_idGlucose, "id_glucose");
ADD_PARAMETER(m_idOxygen, "id_oxygen");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
Anaerobic_Reaction::Anaerobic_Reaction(FEModel* pfem) : FEChemicalReaction(pfem)
{
    m_Rid = m_Pid = -1;
    m_Rtype = false;
    m_region = 0;
}

bool Anaerobic_Reaction::Init()
{
    if (FEChemicalReaction::Init() == false) return false;
    if (m_solR.size() + m_sbmR.size() > 1) {
        feLogError("Provide only one vR for this reaction");
        return false;
    }

    if (m_solP.size() + m_sbmP.size() > 1) {
        feLogError("Provide only one vP for this reaction");
        return false;
    }

    const int ntot = (int)m_v.size();
    for (int itot = 0; itot < ntot; itot++) {
        if (m_vR[itot] > 0) m_Rid = itot;
        if (m_vP[itot] > 0) m_Pid = itot;
    }

    if (m_Rid == -1) {
        feLogError("Provide vR for the reactant");
        return false;
    }
    if (m_Rid >= m_nsol) m_Rtype = true;

    return true;
}

// ============================================================================
// Anaerobic reaction rate
// ============================================================================

double Anaerobic_Reaction::ReactionSupply(FEMaterialPoint& pt)
{
    // get reaction rate
    FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
    double G0 = 1.5; 
    double G = std::min(5.0, std::max(G0, (spt.m_c.empty() ? G0 : spt.m_c[3])));

    double cG = m_psm->GetActualSoluteConcentration(pt, m_idGlucose);
    double cO = m_psm->GetActualSoluteConcentration(pt, m_idOxygen);

    double Ylac = 1.0;
    double alpha_m = 2.0;
    double beta_m = 0.06;
    double delta_m = 0.6;
    double Am = 152.5;
    double cG0 = 25.0;
    double cg_crit = 1.0; //glucose threshold
    double QO2_max = 18.0;
    double KO2 = 0.0127;

    double Qlac_max = 0.0;
    double Ap = 0.0;
    double Bp = 0.0;
    double alpha_p = 0.0;
    double Beta_p = 0.0;

    double rho_cell_eff = 0.0;

    if (m_region == 0) { //NP
        rho_cell_eff = std::min(2000.0, 2000.0 * cG / cg_crit);
        Qlac_max = 160;
        Ap = 1.31;
        Bp = 0.12;
        alpha_p = 35;
        Beta_p = 0.07;
    }
    else if (m_region == 1) { //AF
        rho_cell_eff = std::min(2000.0, 2000.0 * cG / cg_crit);
        Qlac_max = 40;
        Ap = 2.98;
        Bp = 0.28;
        alpha_p = 8;
        Beta_p = 0.02;
    }
    else if (m_region == 2) { //CEP
        rho_cell_eff = std::min(10000.0, 10000.0 * cG / cg_crit);
        Qlac_max = 40;
        Ap = 2.98;
        Bp = 0.28;
        alpha_p = 8;
        Beta_p = 0.02;
    }
    double eta_O2 = delta_m * (1.0 + exp(-cG / alpha_m)) * (1.0 - exp(-cO / beta_m));
    double eta_lac = 1 - eta_O2;
    double Qlac_maint = eta_lac / Ylac * Am * (G - G0);
    double QO2_bio = QO2_max * cO / (KO2 + cO);
    double delta_Qlac = 5 * (QO2_max - QO2_bio);
    double f_nut = 1 - Ap * exp(-cG * cG / alpha_p) * exp(-cO / Beta_p) - Bp * (cG0 - cG) / cG0;
    double Qlac = (Qlac_max * f_nut + delta_Qlac + Qlac_maint) * rho_cell_eff;
    double factor = pow(10, -6) / 3600.0; //unit conversion
    double Qglu = Qlac / 2.0 * factor;
    return Qglu;
}

// ============================================================================
// Anaerobic reaction tangent: deformation
// ============================================================================

mat3ds Anaerobic_Reaction::Tangent_ReactionSupply_Strain(FEMaterialPoint& pt)
{
    double cG = m_psm->GetActualSoluteConcentration(pt, m_idGlucose);
    double cO = m_psm->GetActualSoluteConcentration(pt, m_idOxygen);

    double cGe = m_psm->GetEffectiveSoluteConcentration(pt, m_idGlucose);
    double cOe = m_psm->GetEffectiveSoluteConcentration(pt, m_idOxygen);

    double dcGdJ = m_psm->dkdJ(pt, m_idGlucose) * cGe;
    double dcOdJ = m_psm->dkdJ(pt, m_idOxygen) * cOe;

    FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();

    double G0 = 1.5;
    double G = std::min(5.0, std::max(G0, (spt.m_c.empty() ? G0 : spt.m_c[3])));

    //Parameters
    double QO2_max = 18.0;
    double KO2 = 0.0127;
    double Qlac_max = 0.0;
    double Ap = 0.0;
    double Bp = 0.0;
    double alpha_p = 0.0;
    double Beta_p = 0.0;
    double Ylac = 1.0;
    double alpha_m = 2.0;
    double beta_m = 0.06;
    double delta_m = 0.6;
    double Am = 152.5; 
    double cg_crit = 1.0;
    double cG0 = 25.0;


    double rho_cell_eff = 0.0;
    double drhodcG = 0.0;

    if (m_region == 0) { //NP
        rho_cell_eff = std::min(2000.0, 2000.0 * cG / cg_crit);
        if (cG < cg_crit)
        {
            drhodcG = 2000.0 / cg_crit;
        }
        Qlac_max = 160;
        Ap = 1.31;
        Bp = 0.12;
        alpha_p = 35;
        Beta_p = 0.07;
    }
    else if (m_region == 1) { //AF
        rho_cell_eff = std::min(2000.0, 2000.0 * cG / cg_crit);
        if (cG < cg_crit)
        {
            drhodcG = 2000.0 / cg_crit;
        }
        Qlac_max = 40;
        Ap = 2.98;
        Bp = 0.28;
        alpha_p = 8;
        Beta_p = 0.02;
    }
    else if (m_region == 2) { //CEP
        rho_cell_eff = std::min(10000.0, 10000.0 * cG / cg_crit);
        if (cG < cg_crit)
        {
            drhodcG = 10000.0 / cg_crit;
        }
        Qlac_max = 40;
        Ap = 2.98;
        Bp = 0.28;
        alpha_p = 8;
        Beta_p = 0.02;
    }

    //basal lactate production
    double QO2_bio = QO2_max * cO / (KO2 + cO);
    double delta_Qlac = 5 * (QO2_max - QO2_bio);
    double f_nut = 1 - Ap * exp(-cG * cG / alpha_p) * exp(-cO / Beta_p) - Bp * (cG0 - cG) / cG0;
    double Qlac_bio = (Qlac_max * f_nut + delta_Qlac) * rho_cell_eff;

    //lactate production for cellular maintenance
    double eta_O2 = delta_m * (1.0 + exp(-cG / alpha_m)) * (1.0 - exp(-cO / beta_m));
    double eta_lac = 1 - eta_O2;
    double Qlac_maint = eta_lac / Ylac * Am * (G - G0) * rho_cell_eff;
    
    //derivative basal lactate production
    double dfnutdcO = Ap * exp(-(cG * cG) / alpha_p) * exp(-cO / Beta_p) * (1.0 / Beta_p);
    double dfnutdcG = Ap * exp(-(cG * cG) / alpha_p) * exp(-cO / Beta_p) * (2.0 * cG / alpha_p) + Bp / cG0;
    double dQbiodcO = QO2_max * KO2 / pow(KO2 + cO, 2);
    double ddelta_QlacdcO = -5 * dQbiodcO ;
    double dQlac_biodcO = (Qlac_max * dfnutdcO + ddelta_QlacdcO) * rho_cell_eff;
    double dQlac_biodcG = Qlac_max * dfnutdcG * rho_cell_eff + (Qlac_max * f_nut + delta_Qlac) * drhodcG;

    //derivative lactate production for cellular maintenance
    double dQlac_maintdcO = -delta_m * (1 + exp(-cG / alpha_m)) * 1.0 / beta_m * exp(-cO / beta_m) / Ylac * Am * (G - G0) * rho_cell_eff;
    double dQlac_maintdcG = delta_m * (1 - exp(-cO / beta_m)) * 1.0 / alpha_m * exp(-cG / alpha_m) / Ylac * Am * (G - G0) * rho_cell_eff + eta_lac / Ylac * Am * (G - G0) * drhodcG;
 
    //total derivative
    double dQlacdcG = dQlac_biodcG + dQlac_maintdcG;
    double dQlacdcO = dQlac_biodcO + dQlac_maintdcO;

    double dQlacdJ = dQlacdcG * dcGdJ + dQlacdcO * dcOdJ;
    double factor = pow(10, -6) / 3600.0; //unit conversion
    double dQgludJ = dQlacdJ / 2.0 * factor;

    return mat3dd(1.0) * dQgludJ;
}


// ============================================================================
// Anaerobic reaction tangent: pressure
// ============================================================================

double Anaerobic_Reaction::Tangent_ReactionSupply_Pressure(FEMaterialPoint& pt)
{
    return 0;
}

// ============================================================================
// Anaerobic reaction tangent: effective solute concentration
// ============================================================================

double Anaerobic_Reaction::Tangent_ReactionSupply_Concentration(FEMaterialPoint& pt, int sol)
{
    double cG = m_psm->GetActualSoluteConcentration(pt, m_idGlucose);
    double cO = m_psm->GetActualSoluteConcentration(pt, m_idOxygen);
    double cGe = m_psm->GetEffectiveSoluteConcentration(pt, m_idGlucose);
    double cOe = m_psm->GetEffectiveSoluteConcentration(pt, m_idOxygen);

    double kG = m_psm->GetPartitionCoefficient(pt, m_idGlucose);
    double kO = m_psm->GetPartitionCoefficient(pt, m_idOxygen);

    double dkGdc = m_psm->dkdc(pt, m_idGlucose, m_idGlucose);
    double dkOdc = m_psm->dkdc(pt, m_idOxygen, m_idOxygen);

    FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();

    double G0 = 1.5;
    double G = std::min(5.0, std::max(G0, (spt.m_c.empty() ? G0 : spt.m_c[3])));

    double QO2_max = 18.0;
    double KO2 = 0.0127;

    double Ylac = 1.0;
    double alpha_m = 2.0;
    double beta_m = 0.06;
    double delta_m = 0.6;
    double Am = 152.5;

    double Qlac_max = 0.0;
    double Ap = 0.0;
    double Bp = 0.0;
    double alpha_p = 0.0;
    double Beta_p = 0.0;

    double cg_crit = 1.0;
    double cG0 = 25.0;

    double rho_cell_eff = 0.0;
    double drhodcG = 0.0;

    if (m_region == 0) { //NP
        rho_cell_eff = std::min(2000.0, 2000.0 * cG / cg_crit);
        if (cG < cg_crit)
        {
            drhodcG = 2000.0 / cg_crit;
        }
        Qlac_max = 160;
        Ap = 1.31;
        Bp = 0.12;
        alpha_p = 35;
        Beta_p = 0.07;
    }
    else if (m_region == 1) { //AF
        rho_cell_eff = std::min(2000.0, 2000.0 * cG / cg_crit);
        if (cG < cg_crit)
        {
            drhodcG = 2000.0 / cg_crit;
        }
        Qlac_max = 40;
        Ap = 2.98;
        Bp = 0.28;
        alpha_p = 8;
        Beta_p = 0.02;
    }
    else if (m_region == 2) { //CEP
        rho_cell_eff = std::min(10000.0, 10000.0 * cG / cg_crit);
        if (cG < cg_crit)
        {
            drhodcG = 10000.0 / cg_crit;
        }
        Qlac_max = 40;
        Ap = 2.98;
        Bp = 0.28;
        alpha_p = 8;
        Beta_p = 0.02;
    }
    //basal lactate production
    double QO2_bio = QO2_max * cO / (KO2 + cO);
    double delta_Qlac = 5 * (QO2_max - QO2_bio);
    double f_nut = 1 - Ap * exp(-cG * cG / alpha_p) * exp(-cO / Beta_p) - Bp * (cG0 - cG) / cG0;
    double Qlac_bio = (Qlac_max * f_nut + delta_Qlac) * rho_cell_eff;

    //lactate production for cellular maintenance
    double eta_O2 = delta_m * (1.0 + exp(-cG / alpha_m)) * (1.0 - exp(-cO / beta_m));
    double eta_lac = 1 - eta_O2;
    double Qlac_maint = eta_lac / Ylac * Am * (G - G0) * rho_cell_eff;

    //derivative basal lactate production
    double dfnutdcO = Ap * exp(-(cG * cG) / alpha_p) * exp(-cO / Beta_p) * (1.0 / Beta_p);
    double dfnutdcG = Ap * exp(-(cG * cG) / alpha_p) * exp(-cO / Beta_p) * (2.0 * cG / alpha_p) + Bp / cG0;
    double dQbiodcO = QO2_max * KO2 / pow(KO2 + cO, 2);
    double ddelta_QlacdcO = -5 * dQbiodcO;
    double dQlac_biodcO = (Qlac_max * dfnutdcO + ddelta_QlacdcO) * rho_cell_eff;
    double dQlac_biodcG = Qlac_max * dfnutdcG * rho_cell_eff + (Qlac_max * f_nut + delta_Qlac) * drhodcG;

    //derivative lactate production for cellular maintenance
    double dQlac_maintdcO = -delta_m * (1 + exp(-cG / alpha_m)) * 1.0 / beta_m * exp(-cO / beta_m) / Ylac * Am * (G - G0) * rho_cell_eff;
    double dQlac_maintdcG = delta_m * (1 - exp(-cO / beta_m)) * 1.0 / alpha_m * exp(-cG / alpha_m) / Ylac * Am * (G - G0) * rho_cell_eff + eta_lac / Ylac * Am * (G - G0) * drhodcG;

    //total derivative
    double dQlacdcG = dQlac_biodcG + dQlac_maintdcG;
    double dQlacdcO = dQlac_biodcO + dQlac_maintdcO;

    double dQgludcG = dQlacdcG / 2.0;
    double dQgludcO = dQlacdcO / 2.0;
    double factor = pow(10, -6) / 3600.0; //unit conversion

    if (sol == m_idGlucose)
        return dQgludcG * (kG + cGe * dkGdc) * factor;

    if (sol == m_idOxygen)
        return dQgludcO * (kO + cOe * dkOdc) * factor;

    return 0.0;
}

// ============================================================================
// Degeneration-dependent diffusivity
// ============================================================================

BEGIN_FECORE_CLASS(FEDiffusivity_Calcification, FESoluteDiffusivity)
ADD_PARAMETER(m_free_diff, FE_RANGE_GREATER_OR_EQUAL(0.0), "free_diff")->setUnits(UNIT_DIFFUSIVITY)->setLongName("free diffusivity");
ADD_PARAMETER(m_diff, FE_RANGE_GREATER_OR_EQUAL(0.0), "diff")->setUnits(UNIT_DIFFUSIVITY)->setLongName("diffusivity");
END_FECORE_CLASS();
 
FEDiffusivity_Calcification::FEDiffusivity_Calcification(FEModel* pfem) : FESoluteDiffusivity(pfem)
{
    m_free_diff = 1;
    m_diff = 1;
}

bool FEDiffusivity_Calcification::Validate()
{
    if (FESoluteDiffusivity::Validate() == false) return false;
    FEMesh& mesh = GetMesh();
    for (int i = 0; i < mesh.Elements(); ++i) {
        FEElement& elem = *mesh.Element(i);
        for (int n = 0; n < elem.GaussPoints(); ++n) {
            FEMaterialPoint& mp = *elem.GetMaterialPoint(n);
            if (m_free_diff(mp) < m_diff(mp)) {
                feLogError("free_diff must be >= diff in element %i", elem.GetID());
                return false;
            }
        }
    }
    return true;
}

double FEDiffusivity_Calcification::Free_Diffusivity(FEMaterialPoint& mp)
{
    return m_free_diff(mp);
}

double FEDiffusivity_Calcification::Tangent_Free_Diffusivity_Concentration(FEMaterialPoint& mp, const int isol)
{
    return 0;
}

mat3ds FEDiffusivity_Calcification::Diffusivity(FEMaterialPoint& mp)
{
    FESolutesMaterialPoint& spt = *(mp.ExtractData<FESolutesMaterialPoint>());
    double G = (spt.m_c.empty() ? 0.0 : spt.m_c[3]);
    double calc = (2.0071 * G * G + 6.7071 * G - 3.82) / 100.0; //Equation A.2
    double calc_1_5 = (2.0071 * 1.5 * 1.5 + 6.7071 * 1.5 - 3.82) / 100.0; //initial calcification level (G0=1.5)
    double delta = (1.0 - 1.334 * (calc - calc_1_5)); //Equation A.3
    return mat3dd(delta * m_diff(mp));
}

tens4dmm FEDiffusivity_Calcification::Tangent_Diffusivity_Strain(FEMaterialPoint& mp)
{
    tens4dmm D;
    D.zero();
    return D;
}

mat3ds FEDiffusivity_Calcification::Tangent_Diffusivity_Concentration(FEMaterialPoint& mp, const int isol)
{
    mat3ds d;
    d.zero();
    return d;
}

