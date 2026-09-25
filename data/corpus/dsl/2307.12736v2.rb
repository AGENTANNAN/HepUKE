# BESIII Analysis: e+e- → φη' cross sections at 3.508–4.951 GeV
# Paper: arXiv:2307.12736v2 — two η' decay modes (Rule T1 → separate Algorithm objects)
# Analysis type: Ordinary, multi-energy scan, continuum + search for ψ(3770)→φη'
# 26.1 fb⁻¹ total integrated luminosity across ~41 energy points

### Dataset preparation ###
# Key energy points from Table I across BOSS versions 703, 705, 706, 707
data_3508 = DatasetManager.real_data.find("703_3508")
data_3511 = DatasetManager.real_data.find("703_3511")
data_3582 = DatasetManager.real_data.find("703_3582")
data_3650 = DatasetManager.real_data.find("709_3650")
data_3670 = DatasetManager.real_data.find("703_3670")
data_3773 = DatasetManager.real_data.find("712_3773")
data_3808 = DatasetManager.real_data.find("703_3808")
data_3867 = DatasetManager.real_data.find("703_3867")
data_3871 = DatasetManager.real_data.find("703_3871")
data_3896 = DatasetManager.real_data.find("703_3896")
data_4009 = DatasetManager.real_data.find("703_4009")
data_4085 = DatasetManager.real_data.find("703_4085")
data_4128 = DatasetManager.real_data.find("705_4128")
data_4157 = DatasetManager.real_data.find("705_4157")
data_4178 = DatasetManager.real_data.find("705_4178")
data_4189 = DatasetManager.real_data.find("705_4189")
data_4199 = DatasetManager.real_data.find("705_4199")
data_4209 = DatasetManager.real_data.find("705_4209")
data_4219 = DatasetManager.real_data.find("705_4219")
data_4226 = DatasetManager.real_data.find("705_4226")
data_4236 = DatasetManager.real_data.find("703_4236")
data_4244 = DatasetManager.real_data.find("703_4244")
data_4258 = DatasetManager.real_data.find("703_4258")
data_4267 = DatasetManager.real_data.find("703_4267")
data_4278 = DatasetManager.real_data.find("703_4278")
data_4288 = DatasetManager.real_data.find("703_4288")
data_4308 = DatasetManager.real_data.find("703_4308")
data_4312 = DatasetManager.real_data.find("703_4312")
data_4337 = DatasetManager.real_data.find("703_4337")
data_4358 = DatasetManager.real_data.find("703_4358")
data_4377 = DatasetManager.real_data.find("703_4377")
data_4396 = DatasetManager.real_data.find("703_4396")
data_4416 = DatasetManager.real_data.find("703_4416")
data_4436 = DatasetManager.real_data.find("703_4436")
data_4467 = DatasetManager.real_data.find("703_4467")
data_4527 = DatasetManager.real_data.find("703_4527")
data_4600 = DatasetManager.real_data.find("703_4600")
data_4612 = DatasetManager.real_data.find("706_4612")
data_4628 = DatasetManager.real_data.find("706_4628")
data_4641 = DatasetManager.real_data.find("706_4641")
data_4661 = DatasetManager.real_data.find("706_4661")

# All scan data points
scan_data = [data_3508, data_3511, data_3582, data_3650, data_3670, data_3773,
             data_3808, data_3867, data_3871, data_3896, data_4009, data_4085,
             data_4128, data_4157, data_4178, data_4189, data_4199, data_4209,
             data_4219, data_4226, data_4236, data_4244, data_4258, data_4267,
             data_4278, data_4288, data_4308, data_4312, data_4337, data_4358,
             data_4377, data_4396, data_4416, data_4436, data_4467, data_4527,
             data_4600, data_4612, data_4628, data_4641, data_4661]

# Inclusive MC for key BOSS versions
incMC_703 = DatasetManager.inclusive_mc.find("703_4260")
incMC_705 = DatasetManager.inclusive_mc.find("705_4178")
incMC_706 = DatasetManager.inclusive_mc.find("706_4628")
incMC_709 = DatasetManager.inclusive_mc.find("709_3650")
incMC_712 = DatasetManager.inclusive_mc.find("712_3773")

# ============================================================
# ConExc decay card for e+e- → φη' continuum scan
# ============================================================
# ConExc for continuum production with ISR (Born cross section measurement)
conexc_card = <<~DECAYCARD
    Decay vpho
    1.0000  phi  eta'                     VSS;
    Enddecay

    Decay phi
    1.0000  K+  K-                        VSS;
    Enddecay

    Decay eta'
    1.0000  gamma  pi+  pi-               PHSP;
    Enddecay
End
DECAYCARD

# Exclusive MC for the continuum φη' signal
exMC_signal = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name   = "phi_etap_continuum"
  config.events        = 100_000
  config.decay_card    = conexc_card
  config.cross_section = :default
end

# ============================================================
# Analysis A — Mode I: η' → γπ+π- (4C kinematic fit)
# ============================================================

alg_mode1 = Algorithm.new("PhiEtapMode1_gpipi")
alg_mode1.set_header(["PhiEtapMode1Alg/PhiEtapMode1.h"])

alg_mode1.note(:multi_energy, "Multi-energy scan at 41 energy points from 3.508 to 4.951 GeV; ECMS varies per point")
alg_mode1.note(:conexc_continuum, "Continuum e+e- → φη' production modeled with ConExc; Born cross section extracted via ISR correction")
alg_mode1.note(:vss_model, "φ → K+K- decay modeled with VSS generator")

selection_mode1 = Selection.new
selection_mode1.select_track {
                 cos_theta   0.93
                 Vz   10.0
                 Vr   1.0
                 nTot    "==4"
                 nNet    "==0"
               }
               .pid(method: :probability) {
                 prob_cut   0.001
                 # K/π separation: CL_K > CL_π for kaon candidates
                 identify :kaon, against: [:pion]
                 nkp   "==2"
                 nkm   "==2"
               }
               .remove([:km <= :chrgn])
               .remove([:kp <= :chrgp])
               .assign({:chrgp => :pip, :chrgn => :pim})
               .select_photon {
                 tdc_emc_start   0
                 tdc_emc_end   700
                 angle_to_track   10.0
                 energyThreshold_b 0.025
                 energyThreshold_e  0.050
                 nGam   ">=1"
               }
               .note(:phi_mass_window, "|M(K+K-) - m_φ| < 0.04 GeV/c² applied in ROOT analysis after vertex fit")
               .note(:etap_mass_window, "η' mass window applied in ROOT analysis")
               # 4C kinematic fit for Mode I: constrain total 4-momentum
               .kinematic_fit([:kp, :km, :pip, :pim, :gamma]) {
                 nominal
                 constrain_four_momentum
                 chi2_cut 100
               }

alg_mode1.with_decay_card(conexc_card).apply(selection_mode1)
# NOTE: a competing-hypothesis veto (e.g. against φπ+π-π0) would require
# a second kinematic_fit block without nominal/chi2_cut to store the competing χ² (Rule T2)

# ============================================================
# Analysis B — Mode II: η' → ηπ+π- with η → γγ (5C kinematic fit)
# ============================================================

alg_mode2 = Algorithm.new("PhiEtapMode2_etapipi")
alg_mode2.set_header(["PhiEtapMode2Alg/PhiEtapMode2.h"])

alg_mode2.note(:multi_energy, "Multi-energy scan at 41 energy points from 3.508 to 4.951 GeV; ECMS varies per point")
alg_mode2.note(:conexc_continuum, "Continuum e+e- → φη' production modeled with ConExc; Born cross section extracted via ISR correction")

selection_mode2 = Selection.new
selection_mode2.select_track {
                 cos_theta   0.93
                 Vz   10.0
                 Vr   1.0
                 nTot    "==4"
                 nNet    "==0"
               }
               .pid(method: :probability) {
                 prob_cut   0.001
                 identify :kaon, against: [:pion]
                 nkp   "==2"
                 nkm   "==2"
               }
               .remove([:km <= :chrgn])
               .remove([:kp <= :chrgp])
               .assign({:chrgp => :pip, :chrgn => :pim})
               .select_photon {
                 tdc_emc_start   0
                 tdc_emc_end   700
                 angle_to_track   10.0
                 energyThreshold_b 0.025
                 energyThreshold_e  0.050
                 nGam   ">=2"
               }
               .kalman_kinematic_fit([:gamma, :gamma]) {
                 invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
                 chi2_cut 200
                 neta  ">=1"
               }
               .note(:phi_mass_window, "|M(K+K-) - m_φ| < 0.04 GeV/c² applied in ROOT analysis after vertex fit")
               .note(:etap_mass_window, "η' mass window applied in ROOT analysis")
               # 5C kinematic fit for Mode II: 4C + η mass constraint
               .kinematic_fit([:kp, :km, :pip, :pim, :eta]) {
                 nominal
                 constrain_four_momentum
                 chi2_cut 200
               }

alg_mode2.with_decay_card(conexc_card).apply(selection_mode2)

# ============================================================
# Execute both mode analyses
# ============================================================
alg_mode1.execute_on(scan_data + [incMC_703, incMC_705, incMC_706, incMC_709, incMC_712] + exMC_signal)
alg_mode2.execute_on(scan_data + [incMC_703, incMC_705, incMC_706, incMC_709, incMC_712] + exMC_signal)