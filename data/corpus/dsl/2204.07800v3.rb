### Dataset description ###
# Paper: Cross sections of e+e- -> K+K-J/psi at sqrt(s)=4.127-4.600 GeV
# arXiv: 2204.07800v3
# 28 energy points, 15.6 fb^{-1} total, partial reconstruction (missing K-)
# Two resonant structures observed: Y(4230) and Y(4500)

# 28 energy points from 4.127 to 4.600 GeV
data_4127  = DatasetManager.real_data.find("703_4127")
data_4157  = DatasetManager.real_data.find("703_4157")
data_4178  = DatasetManager.real_data.find("703_4178")
data_4189  = DatasetManager.real_data.find("703_4189")
data_4199  = DatasetManager.real_data.find("703_4199")
data_4209  = DatasetManager.real_data.find("703_4209")
data_4219  = DatasetManager.real_data.find("703_4219")
data_4226  = DatasetManager.real_data.find("703_4226")
data_4236  = DatasetManager.real_data.find("703_4236")
data_4242  = DatasetManager.real_data.find("703_4242")
data_4244  = DatasetManager.real_data.find("703_4244")
data_4258  = DatasetManager.real_data.find("703_4258")
data_4267  = DatasetManager.real_data.find("703_4267")
data_4278  = DatasetManager.real_data.find("703_4278")
data_4287  = DatasetManager.real_data.find("703_4287")
data_4308  = DatasetManager.real_data.find("703_4308")
data_4311  = DatasetManager.real_data.find("703_4311")
data_4337  = DatasetManager.real_data.find("703_4337")
data_4358  = DatasetManager.real_data.find("703_4358")
data_4377  = DatasetManager.real_data.find("703_4377")
data_4387  = DatasetManager.real_data.find("703_4387")
data_4395  = DatasetManager.real_data.find("703_4395")
data_4416  = DatasetManager.real_data.find("703_4416")
data_4436  = DatasetManager.real_data.find("703_4436")
data_4467  = DatasetManager.real_data.find("703_4467")
data_4527  = DatasetManager.real_data.find("703_4527")
data_4574  = DatasetManager.real_data.find("703_4574")
data_4600  = DatasetManager.real_data.find("703_4600")
all_data = [data_4127, data_4157, data_4178, data_4189, data_4199,
            data_4209, data_4219, data_4226, data_4236, data_4242,
            data_4244, data_4258, data_4267, data_4278, data_4287,
            data_4308, data_4311, data_4337, data_4358, data_4377,
            data_4387, data_4395, data_4416, data_4436, data_4467,
            data_4527, data_4574, data_4600]

incMC_4127  = DatasetManager.inclusive_mc.find("703_4127")
incMC_4226  = DatasetManager.inclusive_mc.find("703_4226")
incMC_4436  = DatasetManager.inclusive_mc.find("703_4436")
incMC_4600  = DatasetManager.inclusive_mc.find("703_4600")
all_incMC = [incMC_4127, incMC_4226, incMC_4436, incMC_4600]

# Decay card: e+e- -> K+ K- J/psi, J/psi -> l+l- (two lepton modes)
# Partial reconstruction: only one K+ identified; missing K- inferred via 1C kinematic fit
# Three intermediate components: PHSP (non-resonant), f0(980)J/psi, f2(1270)J/psi
# Electron mode
decay_card_ee = <<~DECAYCARD
    Decay J/psi
    1.0000 e+ e- PHSP;
    Enddecay

    End
DECAYCARD

# Muon mode
decay_card_mumu = <<~DECAYCARD
    Decay J/psi
    1.0000 mu+ mu- PHSP;
    Enddecay

    End
DECAYCARD

# Signal MC: K+K-J/psi with three intermediate components weighted
# Generated at representative energies across the scan
exMC_ee = DatasetManager.create_exclusive_mc_for([data_4226, data_4436, data_4600]) do |config|
  config.sample_name = "exmc_kkjpsi_ee"
  config.events = 300000
  config.decay_card = decay_card_ee
  config.cross_section = :default
end

exMC_mumu = DatasetManager.create_exclusive_mc_for([data_4226, data_4436, data_4600]) do |config|
  config.sample_name = "exmc_kkjpsi_mumu"
  config.events = 300000
  config.decay_card = decay_card_mumu
  config.cross_section = :default
end

### Event selection (BOSS) — J/psi -> e+e- mode ###
# Partial reconstruction: reconstruct e+e- (J/psi) and at least one K+
# Missing K- inferred via 1C kinematic fit on kaon mass
alg_name_ee = "KKJpsiEE"
alg_ee = Algorithm.new(alg_name_ee)
alg_ee.set_header(["#{alg_name_ee}Alg/#{alg_name_ee}.h"])
      .set_constant({ "ECMS" => [:double, 4.226] })

selection_ee = Selection.new
selection_ee.select_track {
               cos_theta 0.93
               Vz 10.0
               Vr 1.0
               nTot ">=3"       # e+, e-, and at least one K+
             }
             .pid(method: :probability) {
               identify :ep, :em, against: :kp
               identify_high_momentum_leptons :ep, :em
               prob_cut 0.001
             }
             .note(:lepton_identification, "Tracks with p > 1.0 GeV assigned as leptons; E/p > 0.8 for e+/e-")
             .pid(method: :probability) {
               identify :kp, against: :pip
               prob_cut 0.001
             }
             .note(:kaon_pid, "C.L.(K) > 0.001 and C.L.(K) > C.L.(pi) for each kaon candidate")
             .kinematic_fit([:ep, :em, :kp]) {
               nominal
               constrain_four_momentum
               chi2_cut 200
             }
             .note(:partial_reconstruction, "Partial reconstruction: only one K+ identified; missing K- inferred via 1C kinematic fit constraining missing mass to known kaon mass; chi2 = vertex_fit_chi2 + 1C_fit_chi2 < 20; if multiple K+ candidates, choose minimum chi2")
             .note(:bhabha_veto, "Radiative Bhabha background removal: all opposite-charge track pairs required to have cos(theta_open) < 0.98")
             .note(:intermediate_states, "Three components weighted: PHSP (non-resonant K+K-J/psi), f0(980)J/psi, f2(1270)J/psi; weights extracted from fit to M(K+K-) distribution at high-yield energies")
             .note(:jpsi_signal, "J/psi signal region: M(l+l-) within 2.5*resolution; sidebands: same-width intervals shifted by 0.01 GeV from signal boundaries")
             .note(:born_cross_section, "Born cross section extracted in ROOT: sigma_B = N_obs / (L_int * epsilon * (1+delta)_ISR * 1/|1-Pi|^2 * B(J/psi->l+l-)); ISR correction factor from iterative method")

alg_ee.with_decay_card(decay_card_ee).apply(selection_ee)
root_files_ee = alg_ee.execute_on(all_data + all_incMC + exMC_ee)

### Event selection (BOSS) — J/psi -> mu+mu- mode ###
alg_name_mumu = "KKJpsiMuMu"
alg_mumu = Algorithm.new(alg_name_mumu)
alg_mumu.set_header(["#{alg_name_mumu}Alg/#{alg_name_mumu}.h"])
        .set_constant({ "ECMS" => [:double, 4.226] })

selection_mumu = Selection.new
selection_mumu.select_track {
                 cos_theta 0.93
                 Vz 10.0
                 Vr 1.0
                 nTot ">=3"       # mu+, mu-, and at least one K+
               }
               .pid(method: :probability) {
                 identify :mup, :mum, against: :kp
                 identify_high_momentum_leptons :mup, :mum
                 prob_cut 0.001
               }
               .note(:muon_identification, "Tracks with p > 1.0 GeV assigned as leptons; E/p < 0.8 for mu+/mu-")
               .pid(method: :probability) {
                 identify :kp, against: :pip
                 prob_cut 0.001
               }
               .note(:kaon_pid, "C.L.(K) > 0.001 and C.L.(K) > C.L.(pi) for each kaon candidate")
               .kinematic_fit([:mup, :mum, :kp]) {
                 nominal
                 constrain_four_momentum
                 chi2_cut 200
               }
               .note(:partial_reconstruction, "Partial reconstruction: only one K+ identified; missing K- inferred via 1C kinematic fit; chi2 < 20; best K+ candidate chosen by minimum total chi2")
               .note(:muon_counter, "Muon counter penetration depth > 40 cm for at least one muon candidate to suppress hadron backgrounds")
               .note(:intermediate_states, "Three components weighted: PHSP, f0(980)J/psi, f2(1270)J/psi")
               .note(:jpsi_signal, "J/psi signal region: M(l+l-) within 2.5*resolution; sidebands: same-width intervals shifted by 0.01 GeV")
               .note(:born_cross_section, "Born cross section extracted in ROOT: sigma_B = N_obs / (L_int * epsilon * (1+delta)_ISR * 1/|1-Pi|^2 * B(J/psi->l+l-))")

alg_mumu.with_decay_card(decay_card_mumu).apply(selection_mumu)
root_files_mumu = alg_mumu.execute_on(all_data + all_incMC + exMC_mumu)