# 2206.08554v2: Resonance structures in e+e- -> pi+ pi- J/psi
#
# Ordinary analysis. Cross section measurement at 40+ energy points
# from 3.7730 to 4.7008 GeV (~23 fb-1 total).
# Two independent J/psi decay modes: J/psi -> e+e- and J/psi -> mu+mu-.
# Event selection: 4 charged tracks, net charge zero.
# Lepton ID: high-momentum (p > 1 GeV) tracks with e/mu separation via EMC.
# 4C kinematic fit with chi2 < 60.
# Pion dE/dx discriminator. BDT for two-photon background (e mode).
# cos opening angle cuts. ISR correction iterative procedure (ROOT side).

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

# --- Datasets: XYZ data (40+ points) + R-scan data (13 points) ---
# Main energy points used in the cross section scan
data_3773  = DatasetManager.real_data.find("712_3773")
data_4130  = DatasetManager.real_data.find("705_4130")
data_4160  = DatasetManager.real_data.find("705_4160")
data_4180  = DatasetManager.real_data.find("703_4180")
data_4190  = DatasetManager.real_data.find("703_4190")
data_4200  = DatasetManager.real_data.find("703_4200")
data_4210  = DatasetManager.real_data.find("703_4210")
data_4220  = DatasetManager.real_data.find("703_4220")
data_4230  = DatasetManager.real_data.find("703_4230")
data_4237  = DatasetManager.real_data.find("703_4237")
data_4246  = DatasetManager.real_data.find("703_4246")
data_4260  = DatasetManager.real_data.find("703_4260")
data_4270  = DatasetManager.real_data.find("703_4270")
data_4280  = DatasetManager.real_data.find("703_4280")
data_4290  = DatasetManager.real_data.find("705_4290")
data_4315  = DatasetManager.real_data.find("705_4315")
data_4340  = DatasetManager.real_data.find("705_4340")
data_4360  = DatasetManager.real_data.find("703_4360")
data_4380  = DatasetManager.real_data.find("705_4380")
data_4400  = DatasetManager.real_data.find("705_4400")
data_4420  = DatasetManager.real_data.find("703_4420")
data_4440  = DatasetManager.real_data.find("705_4440")
data_4470  = DatasetManager.real_data.find("703_4470")
data_4530  = DatasetManager.real_data.find("703_4530")
data_4575  = DatasetManager.real_data.find("703_4575")
data_4600  = DatasetManager.real_data.find("703_4600")
data_4610  = DatasetManager.real_data.find("706_4610")
data_4620  = DatasetManager.real_data.find("706_4620")
data_4640  = DatasetManager.real_data.find("706_4640")
data_4660  = DatasetManager.real_data.find("706_4660")
data_4680  = DatasetManager.real_data.find("706_4680")
data_4700  = DatasetManager.real_data.find("706_4700")

# R-scan data (additional scan points 4.41-4.59 GeV)
data_rscan_4420 = DatasetManager.real_data.find("703_4420") rescue nil
data_rscan_4470 = DatasetManager.real_data.find("703_4470") rescue nil
data_rscan_4530 = DatasetManager.real_data.find("703_4530") rescue nil
data_rscan_4575 = DatasetManager.real_data.find("703_4575") rescue nil

xyz_data = [data_3773, data_4130, data_4160, data_4180, data_4190,
            data_4200, data_4210, data_4220, data_4230, data_4237,
            data_4246, data_4260, data_4270, data_4280, data_4290,
            data_4315, data_4340, data_4360, data_4380, data_4400,
            data_4420, data_4440, data_4470, data_4530, data_4575,
            data_4600, data_4610, data_4620, data_4640, data_4660,
            data_4680, data_4700]

# ============================================
# Decay cards (one per J/psi decay mode)
# ============================================

decay_card_ee = <<~DECAYCARD
  Decay psi(4260)
  1.0 pi+ pi- J/psi PHSP;
  Enddecay
  Decay J/psi
  1.0 e+ e- PHOTOS VLL;
  Enddecay
  End
DECAYCARD

decay_card_mumu = <<~DECAYCARD
  Decay psi(4260)
  1.0 pi+ pi- J/psi PHSP;
  Enddecay
  Decay J/psi
  1.0 mu+ mu- PHOTOS VLL;
  Enddecay
  End
DECAYCARD

# ============================================
# Algorithm 1: J/psi -> e+ e- mode
# ============================================

alg_ee = Algorithm.new("PiPiJPsiEE")
alg_ee.set_header(["PiPiJPsiEEAlg/PiPiJPsiEE.h"])

# Common event selection for e+e- mode
event_selection_ee = Selection.new
  .select_track do
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nTot      "==4"
    nNet      "==0"
  end
  .pid(method: :probability) do
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6
    identify :pion, against: [:kaon, :proton]
    nlp "==1"
    nlm "==1"
    npip "==1"
    npim "==1"
  end
  .kinematic_fit([:pip, :pim, :lp, :lm]) do
    constrain_four_momentum
    chi2_cut 60
    nominal
  end

alg_ee.note(:cos_opening_angle,
  "cos(pi+,pi-) < 0.98 and cos(pi+-,e-+) < 0.98 to suppress " \
  "gamma conversion from radiative Bhabha and di-muon backgrounds. " \
  "Applied via for_each filter after kinematic fit in ROOT stage.")

alg_ee.note(:pion_dedx_discriminator,
  "Pion PID via dE/dx discriminator: chi_pi+ < 2.8, chi_pi- < 3.0 " \
  "for the e+e- mode. Reduces low-momentum electron background. " \
  "Applied as additional pid cut in ROOT stage.")

alg_ee.note(:bdt_background_suppression,
  "BDT trained using e+e- -> e+e- mu+mu- (two-photon process) MC " \
  "to suppress two-photon background. Input variables: EMC energy, " \
  "TOF, dE/dx, opening angles. BDT response cut optimized via " \
  "S/sqrt(S+B). Applied in ROOT analysis stage.")

alg_ee.note(:isr_correction,
  "ISR correction factor (1+delta) computed via iterative procedure " \
  "using KKMC. Flat cross section assumed initially, iterated until " \
  "convergence < 0.1%. XYZ data used for iteration; R-scan factors " \
  "derived from the fitted model. Applied in ROOT stage.")

alg_ee.note(:em_shower_separation,
  "Electron-muon separation in EMC: electron candidates require " \
  "EMC energy > 1.1 GeV; muon candidates require EMC energy < 0.4 GeV. " \
  "Applied in the lepton identification step. Further refined in ROOT.")

alg_ee.with_decay_card(decay_card_ee).apply(event_selection_ee)
alg_ee.execute_on(xyz_data)

# ============================================
# Algorithm 2: J/psi -> mu+ mu- mode
# ============================================

alg_mumu = Algorithm.new("PiPiJPsiMuMu")
alg_mumu.set_header(["PiPiJPsiMuMuAlg/PiPiJPsiMuMu.h"])

event_selection_mumu = Selection.new
  .select_track do
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nTot      "==4"
    nNet      "==0"
  end
  .pid(method: :probability) do
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6
    identify :pion, against: [:kaon, :proton]
    nlp "==1"
    nlm "==1"
    npip "==1"
    npim "==1"
  end
  .kinematic_fit([:pip, :pim, :lp, :lm]) do
    constrain_four_momentum
    chi2_cut 60
    nominal
  end

alg_mumu.note(:cos_opening_angle,
  "cos(pi+,pi-) < 0.98 to suppress gamma conversion background. " \
  "The pion-electron opening angle cut is not applied for the mu+mu- mode. " \
  "Applied via for_each filter after kinematic fit in ROOT stage.")

alg_mumu.note(:pion_dedx_discriminator,
  "Pion PID via dE/dx discriminator: chi_pi+ < 3.2, chi_pi- < 4.10 " \
  "for the mu+mu- mode. Applied as additional pid cut in ROOT stage.")

alg_mumu.note(:em_shower_separation,
  "Electron-muon separation in EMC: muon candidates require " \
  "EMC energy < 0.4 GeV; electron candidates require EMC energy > 1.1 GeV. " \
  "Applied in the lepton identification step.")

alg_mumu.note(:isr_correction,
  "ISR correction factor (1+delta) computed via same iterative procedure " \
  "as the e+e- mode. Cross section given by weighted average of e+e- and " \
  "mu+mu- modes (inverse statistical uncertainty weights). " \
  "Applied in ROOT stage.")

alg_mumu.note(:cross_section_averaging,
  "Final cross section = weighted average of e+e- and mu+mu- mode results, " \
  "using inverse statistical uncertainties as weights. " \
  "Performed in ROOT analysis. Both modes combined after independent " \
  "signal yield extraction from M(l+l-) fits.")

alg_mumu.with_decay_card(decay_card_mumu).apply(event_selection_mumu)
alg_mumu.execute_on(xyz_data)