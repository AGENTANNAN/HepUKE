# BESIII analysis: Measurement of e+e- -> K+ K- psi(2S) cross section at sqrt(s)=4.699-4.951 GeV
# arXiv: 2407.20009v1
# Partial reconstruction with 4 approaches, search for Z_cs+- in K+- psi(2S)
# 2.5 fb^-1 at 6 energy points: 4.700, 4.720, 4.740, 4.750, 4.780, 4.843, 4.914, 4.946

### Datasets ###
data_4700 = DatasetManager.real_data.find("706_4700")
data_4720 = DatasetManager.real_data.find("707_4720")
data_4740 = DatasetManager.real_data.find("707_4740")
data_4750 = DatasetManager.real_data.find("707_4750")
data_4780 = DatasetManager.real_data.find("707_4780")
data_4840 = DatasetManager.real_data.find("707_4840")
data_4914 = DatasetManager.real_data.find("707_4914")
data_4946 = DatasetManager.real_data.find("707_4946")

all_energy_points = [data_4700, data_4720, data_4740, data_4750,
                     data_4780, data_4840, data_4914, data_4946]

incMC_4700 = DatasetManager.inclusive_mc.find("706_4700")
incMC_4720 = DatasetManager.inclusive_mc.find("707_4720")
incMC_4740 = DatasetManager.inclusive_mc.find("707_4740")
incMC_4750 = DatasetManager.inclusive_mc.find("707_4750")
incMC_4780 = DatasetManager.inclusive_mc.find("707_4780")
incMC_4840 = DatasetManager.inclusive_mc.find("707_4840")
incMC_4914 = DatasetManager.inclusive_mc.find("707_4914")
incMC_4946 = DatasetManager.inclusive_mc.find("707_4946")

all_incMC = [incMC_4700, incMC_4720, incMC_4740, incMC_4750,
             incMC_4780, incMC_4840, incMC_4914, incMC_4946]

### Decay cards ###
# Approach (i): tag K+ K- J/psi, partial reco psi(2S) via RM(K+K-)
# psi(2S) -> J/psi + X, J/psi -> l+ l- (l = e, mu)
decay_card_approach_i = <<~DECAYCARD
    Decay vpho
    1.000 K+ K- psi(2S) PHSP;
    Enddecay

    Decay psi(2S)
    1.000 J/psi pi+ pi- PHSP;
    Enddecay

    Decay J/psi
    1.000 e+ e- PHSP;
    Enddecay

    End
DECAYCARD

# Approach (ii): tag one K+-, reconstruct psi(2S) -> pi+ pi- J/psi, 1C fit to missing kaon
decay_card_approach_ii = <<~DECAYCARD
    Decay vpho
    1.000 K+ K- psi(2S) PHSP;
    Enddecay

    Decay psi(2S)
    1.000 J/psi pi+ pi- PHSP;
    Enddecay

    Decay J/psi
    1.000 mu+ mu- PHSP;
    Enddecay

    End
DECAYCARD

# Approach (iii): tag K+ K- psi(2S) with psi(2S) -> l+ l-
decay_card_approach_iii = <<~DECAYCARD
    Decay vpho
    1.000 K+ K- psi(2S) PHSP;
    Enddecay

    Decay psi(2S)
    1.000 e+ e- PHSP;
    Enddecay

    End
DECAYCARD

# Approach (iv): tag one K+-, reconstruct psi(2S) -> l+ l-, 1C fit to missing kaon
decay_card_approach_iv = <<~DECAYCARD
    Decay vpho
    1.000 K+ K- psi(2S) PHSP;
    Enddecay

    Decay psi(2S)
    1.000 mu+ mu- PHSP;
    Enddecay

    End
DECAYCARD

### Exclusive MC ###
exMC_all = DatasetManager.create_exclusive_mc_for(all_energy_points) do |config|
  config.sample_name   = "ee_kkpsip_signal"
  config.events        = 500_000
  config.decay_card    = decay_card_approach_i
  config.cross_section = :default
end

### Algorithm: e+e- -> K+ K- psi(2S) cross section (4 orthogonal partial-reco approaches) ###
alg_cross_section = Algorithm.new("EEK KPsi2SCrossSection")
alg_cross_section.set_header(["EEK KPsi2SAlg/EEK KPsi2S.h"])
  .set_constant({"ECMS" => [:double, 4.843]})
  .note(:four_approaches, "Four orthogonal partial reconstruction approaches:
    (i) tag both K+, K-, and J/psi from psi(2S)->J/psi+X; signal via RM(K+K-).
    (ii) tag one K+-, psi(2S)->pi+pi-J/psi; 1C fit to missing kaon mass (chi2_1C < 50).
    (iii) tag both K+, K-, psi(2S)->l+l-.
    (iv) tag one K+-, psi(2S)->l+l-; 1C fit to missing kaon mass (chi2_1C < 15).
    Approaches are orthogonal by construction.")
  .note(:charged_track, "Same criteria as Ref. [14]: |cos_theta|<0.93, V_xy<1 cm, |V_z|<10 cm. Kaon PID: CL_K > CL_pi; Pion PID: CL_pi > CL_K.")
  .note(:lepton_id, "Electrons: E/p > 0.8. Muons: E/p < 0.4. J/psi->mumu requires >=1 muon penetrating >3 MUC layers.")
  .note(:lepton_pair, "Two charged tracks with p > 1.0 GeV/c and opposite charge as lepton pair from J/psi or psi(2S).")
  .note(:jpsi_window, "J/psi mass window: M(l+l-) in (3.05, 3.15) GeV/c^2.")
  .note(:psip_window, "psi(2S) mass window for approach (iii): M(l+l-) in (3.631, 3.726) GeV/c^2.")
  .note(:bhabha_veto, "Approach (iv) psi(2S)->e+e- channel: veto cos(theta_e+) > 0.85 and cos(theta_e-) < -0.85.")
  .note(:rm_kk_signal, "psi(2S) signal extraction from RM(K+K-) distribution. Signal region [3.67, 3.71] GeV/c^2. Sidebands 2x wider.")
  .note(:signal_counting, "Signal yield from Poisson likelihood: L(x,y|s,b,tau) = Pois(x|s+tau*b)*Pois(y|b). Profile likelihood for errors and significance.")
  .note(:cross_section, "Born cross section calculated with ISR correction via iterative BW lineshape, vacuum polarization factor 1.055.")
  .note(:zcs_search, "Z_cs+- search via simultaneous fit to RM^2(K+) and RM^2(K-) at 4.843 GeV. Z_cs modeled with BW*PS integral convolved with resolution.")
  .note(:partial_reco_missing_k, "Approaches (ii) and (iv): one kaon missing, reconstructed via 1C kinematic fit constraining missing kaon mass. This is handled in the BOSS analysis framework.")

sel_cross_section = Selection.new
sel_cross_section.select_track {
                     cos_theta 0.93
                     Vz 10.0
                     Vr 1.0
                     nChrp ">=1"
                     nChrn ">=1"
                   }
                  .select_photon {
                     tdc_emc_start 0
                     tdc_emc_end 700
                     angle_to_track 10.0
                     energyThreshold_b 0.025
                     energyThreshold_e 0.050
                     nGam ">=0"
                   }
                  .pid(method: :probability) {
                     prob_cut 0.001
                     identify :kaon, against: [:pion]
                     identify :pion, against: [:kaon]
                     nkp ">=0"
                     nkm ">=0"
                     npip ">=0"
                     npim ">=0"
                   }
                  .kinematic_fit([:kp, :km, :ep, :em, :pip, :pim]) {
                     nominal
                     constrain_four_momentum
                     chi2_cut 200
                   }

alg_cross_section.with_decay_card(decay_card_approach_i).apply(sel_cross_section)
alg_cross_section.execute_on(all_energy_points + all_incMC + exMC_all)