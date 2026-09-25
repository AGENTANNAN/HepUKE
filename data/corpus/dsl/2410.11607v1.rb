# BESIII observation of χ_cJ → p pbar K_S^0 K^- π^+ + c.c. (J=0,1,2)
# ArXiv: 2410.11607v1
# Dataset: (27.12±0.14)×10^8 ψ(3686) events

### Dataset preparation ###
psip_data = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# Decay card: ψ(3686) → γ χ_c0; χ_c0 → p+ anti-p- K_S0 K- π+
# (representative; χ_c1 and χ_c2 use analogous cards)
decay_card = <<~DECAYCARD
    Decay psi(2S)
    1.0000  gamma  chi_c0  PHSP;
    Enddecay

    Decay chi_c0
    1.0000  p+  anti-p-  K_S0  K-  pi+  PHSP;
    Enddecay

    Decay K_S0
    1.0000  pi+  pi-  PHSP;
    Enddecay

    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "psip_gamma_chi_c0_ppbar_Ks_K_pi"
  config.related_dataset = psip_data
  config.events = 1_000_000
  config.decay_card = decay_card
  config.cross_section = :default
end

### Event selection: χ_cJ → p pbar K_S^0 K^- π^+ + c.c. ###
alg = Algorithm.new("ChiCJToPPbarKsKPi")
alg.set_header(["ChiCJToPPbarKsKPiAlg/ChiCJToPPbarKsKPi.h"])
  .set_constant(ECMS: 3.686)

event_selection = Selection.new
  .select_track {
    cos_theta  0.93
    Vz         10.0
    Vr         1.0
    nChrp      ">=3"
    nChrn      ">=3"
    nNet       "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
    nGam              ">=1"
  }
  .pid(method: :probability) {
    prob_cut   0.001
    identify :proton, against: [:kaon, :pion]
    nprp  ">=1"
    nprm  ">=1"
    identify :kaon, against: [:pion]
  }
  .remove([:prp <= :chrgp])
  .remove([:prm <= :chrgn])
  .remove([:km <= :chrgn])
  .remove([:kp <= :chrgp])
  .assign({:chrgp => :pip, :chrgn => :pim})
  .secondary_vertex_fit([:pip, :pim]) do
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  # Main mode: χ_cJ → p pbar K_S^0 K^- π^+ (nominal)
  .kinematic_fit([:prp, :prm, :km, :pip, :K_S0, :gamma]) {
    nominal
    constrain_four_momentum
    chi2_cut 50
  }
  # Charge-conjugate mode: χ_cJ → anti-p p K_S^0 K^+ π^-
  .kinematic_fit([:prp, :prm, :kp, :pim, :K_S0, :gamma]) {
    constrain_four_momentum
    chi2_cut 50
  }

alg.note(:ks_reconstruction, "K_S0 from two oppositely charged tracks; first vertex fit chi2<200; second vertex fit requires decay length significance >2; mass window |M(pi+pi-)-m(K_S0)|<0.012 GeV/c^2 applied in ROOT analysis")
  .note(:ks_sideband, "K_S0 sideband region 0.020 < |M(pi+pi-)-m(K_S0)| < 0.044 GeV/c^2 used for background subtraction in the 2D simultaneous fit")
  .note(:chi_cJ_states, "three chi_cJ states (J=0,1,2) analyzed; this DSL uses chi_c0 as representative in the decay card; chi_c1 and chi_c2 use analogous cards")
  .note(:mixed_signal_mc, "signal MC is a mixture of: (a) PHSP chi_cJ→ppbar K_S0 K pi, (b) chi_cJ→pbar Lambda(1520) K_S0 pi with Lambda(1520)→p K-, and (c) chi_cJ→ppbar K*+ K- with K*+→K_S0 pi+; mixing fractions from data 2D fit")
  .note(:chi2_optimization, "chi2_4C<50 optimized with FOM=S/sqrt(S+B) using signal MC and inclusive MC")
  .note(:pwa_fit, "branching fractions extracted via 1D simultaneous unbinned ML fit to M(ppbar K_S0 K pi) in K_S0 signal and sideband regions; Lambda(1520) yields from 2D simultaneous fit to M(ppbar K_S0 K pi) vs M(pK)")
  .with_decay_card(decay_card)
  .apply(event_selection)

alg.execute_on([psip_data, psip_incMC, exMC_signal])