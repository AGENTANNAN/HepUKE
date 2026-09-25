### Dataset description ###
jpsi_data   = DatasetManager.real_data.find("708_3097")     # J/psi real data (225.2e6 events)
jpsi_incMC  = DatasetManager.inclusive_mc.find("708_3097")  # J/psi inclusive MC (Lund model, 2e8 events)
psip_data   = DatasetManager.real_data.find("709_3686")     # psi(2S) real data (1.06e8 events)
psip_incMC  = DatasetManager.inclusive_mc.find("709_3686")  # psi(2S) inclusive MC

### Decay cards ###
decay_card_jpsi_gpp = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma p+ anti-p-                                PHSP;
    Enddecay

    End
DECAYCARD

decay_card_psip_gpp = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma p+ anti-p-                                PHSP;
    Enddecay

    End
DECAYCARD

# Dominant BOSS-side background for J/psi -> gamma p pbar: J/psi -> pi0 p pbar (asymmetric pi0 decays)
decay_card_jpsi_pi0pp = <<~DECAYCARD
    Decay J/psi
    1.0000 pi0 p+ anti-p-                                  PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma                                     PHSP;
    Enddecay

    End
DECAYCARD

decay_card_psip_pi0pp = <<~DECAYCARD
    Decay psi(2S)
    1.0000 pi0 p+ anti-p-                                  PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma                                     PHSP;
    Enddecay

    End
DECAYCARD

### Exclusive MC samples ###
exMC_jpsi_gpp   = DatasetManager.create_exclusive_mc { |c| c.sample_name="jpsi_gamma_ppbar";      c.related_dataset=jpsi_data; c.events=200000; c.decay_card=decay_card_jpsi_gpp;   c.cross_section=:default }
exMC_jpsi_pi0pp = DatasetManager.create_exclusive_mc { |c| c.sample_name="jpsi_pi0_ppbar";        c.related_dataset=jpsi_data; c.events=200000; c.decay_card=decay_card_jpsi_pi0pp; c.cross_section=:default }
exMC_psip_gpp   = DatasetManager.create_exclusive_mc { |c| c.sample_name="psip_gamma_ppbar";      c.related_dataset=psip_data; c.events=200000; c.decay_card=decay_card_psip_gpp;   c.cross_section=:default }
exMC_psip_pi0pp = DatasetManager.create_exclusive_mc { |c| c.sample_name="psip_pi0_ppbar";        c.related_dataset=psip_data; c.events=200000; c.decay_card=decay_card_psip_pi0pp; c.cross_section=:default }


### Event selection — Mode I: J/psi -> gamma p pbar ###
alg_jpsi = Algorithm.new("JpsiGammaPPbar")
alg_jpsi.set_header(["JpsiGammaPPbarAlg/JpsiGammaPPbar.h"])
        .set_constant({"ECMS" => [:double, 3.097]})

sel_jpsi = Selection.new
sel_jpsi.select_track {
           nChrp   "==1"     # one proton
           nChrn   "==1"     # one anti-proton
           nNet    "==0"
           cos_theta 0.93
           Vz      10.0
           Vr      1.0
         }
         .select_photon {
           nGam                ">=1"
           energyThreshold_b   0.025      # >25 MeV barrel (|cos theta|<0.8)
           energyThreshold_e   0.050      # >50 MeV endcap (0.86<|cos theta|<0.92)
           tdc_emc_start       0
           tdc_emc_end         14
         }
         .pid(method: :probability) {
           prob_cut  0.001
           identify :proton, against: [:kaon, :pion]   # assign highest-CL hypothesis (p/anti-p)
         }
         .select_isolated_photon {
           angle_to_prm_track  30.0    # gamma isolated from anti-proton by >30 deg
           angle_to_prp_track  30.0
           nGam                ">=1"
         }
         .kinematic_fit([:gamma, :prp, :prm]) {
           nominal
           constrain_four_momentum
           chi2_cut 20              # chi^2_4C < 20
         }

alg_jpsi.note(:umiss_cut,
              "|U_miss| < 0.05 GeV, where U_miss = E_miss - |P_miss| computed from the two charged tracks; suppresses multi-photon backgrounds.")
        .note(:ptgamma_cut,
              "P_t_gamma^2 < 0.0005 (GeV/c)^2 with P_t_gamma^2 = 4|P_miss|^2 sin^2(theta_gamma/2); further suppresses multi-photon backgrounds.")
        .note(:low_momentum_track_veto,
              "reject events containing any track with momentum below 0.3 GeV/c (data/MC efficiency discrepancy for soft tracks).")

alg_jpsi.with_decay_card(decay_card_jpsi_gpp).apply(sel_jpsi)
alg_jpsi.execute_on([jpsi_data, jpsi_incMC, exMC_jpsi_gpp, exMC_jpsi_pi0pp])


### Event selection — Mode II: psi(2S) -> gamma p pbar (analogous selection) ###
alg_psip = Algorithm.new("PsipGammaPPbar")
alg_psip.set_header(["PsipGammaPPbarAlg/PsipGammaPPbar.h"])
        .set_constant({"ECMS" => [:double, 3.686]})

sel_psip = Selection.new
sel_psip.select_track {
           nChrp   "==1"
           nChrn   "==1"
           nNet    "==0"
           cos_theta 0.93
           Vz      10.0
           Vr      1.0
         }
         .select_photon {
           nGam                ">=1"
           energyThreshold_b   0.025
           energyThreshold_e   0.050
           tdc_emc_start       0
           tdc_emc_end         14
         }
         .pid(method: :probability) {
           prob_cut  0.001
           identify :proton, against: [:kaon, :pion]
         }
         .select_isolated_photon {
           angle_to_prm_track  30.0
           angle_to_prp_track  30.0
           nGam                ">=1"
         }
         .kinematic_fit([:gamma, :prp, :prm]) {
           nominal
           constrain_four_momentum
           chi2_cut 20
         }

alg_psip.note(:umiss_cut,
              "|U_miss| < 0.05 GeV; same definition as in J/psi mode.")
        .note(:ptgamma_cut,
              "P_t_gamma^2 < 0.0005 (GeV/c)^2; same definition as in J/psi mode.")
        .note(:low_momentum_track_veto,
              "reject events containing any track with momentum below 0.3 GeV/c.")

alg_psip.with_decay_card(decay_card_psip_gpp).apply(sel_psip)
alg_psip.execute_on([psip_data, psip_incMC, exMC_psip_gpp, exMC_psip_pi0pp])
