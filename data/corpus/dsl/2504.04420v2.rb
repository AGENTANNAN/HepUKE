# Paper: 2504.04420v2 — First observation of psi(3686) -> Xi- K_S0 Omega+ + c.c.
# Partial reconstruction: Xi- and K_S0 reconstructed, Omega+ via recoil mass RM(Xi- K_S0)
# Sample: psi(3686) at 3.686 GeV

### Dataset preparation ###
psip_data = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

decay_card_signal = <<~DECAYCARD
    Decay psi(3686)
    1.000 Xi- K_S0 Omega+ PHSP;
    Enddecay

    Decay Xi-
    1.000 Lambda pi- PHSP;
    Enddecay

    Decay Lambda
    1.000 p+ pi- HypWK;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_psip_XiKsOmega"
  config.related_dataset = psip_data
  config.events = 2_000_000
  config.decay_card = decay_card_signal
  config.cross_section = :default
end

### Event selection ###
alg = Algorithm.new("XiKsOmega")
alg.set_header(["XiKsOmegaAlg/XiKsOmega.h"])
    .set_constant({"ECMS" => [:double, 3.686]})

event_selection = Selection.new
event_selection
  .select_track {
    cos_theta 0.93
    Vz 20.0        # |Vz| < 20 cm
    Vr 1.0
    nChrp ">=2"    # p+ (Lambda) + pi+ (K_S0) + Omega+ inclusive tracks
    nChrn ">=3"    # pi- (Lambda) + pi- (Xi-) + pi- (K_S0) + Omega+ inclusive tracks
  }
  # PID: identify proton, remove from chrgp, remaining are pions
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprp ">=1"
  }
  .remove([:prp <= :chrgp])
  # Remaining charged tracks assigned as pions
  .assign({:chrgp => :pip, :chrgn => :pim})
  # Reconstruct Lambda -> p pi- via secondary vertex fit
  .secondary_vertex_fit([:prp, :pim]) {
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  # Reconstruct K_S0 -> pi+ pi- via secondary vertex fit
  .secondary_vertex_fit([:pip, :pim]) {
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  # Partial reconstruction: miss Omega+, reconstruct Xi- + K_S0
  .partial_miss([3]) {
    # Omega+ peak fitted in ROOT; no require_recoil_mass cut applied
    # Xi- mass window [1.313, 1.330] applied in ROOT
  }

# Xi- reconstruction specifics (vertex fit, decay length, best candidate) not expressible in DSL
alg.note(:Xi_reconstruction,
  "Xi- reconstructed via vertex fit on Lambda + pi-; decay length d_Xi > 0 required; best candidate by minimum sqrt((M_Lambda_pi-M_Xi_PDG)/sigma_M)^2+((M_p_pi-M_Lambda_PDG)/sigma_M)^2; mass window M_Xi in [1.313, 1.330] GeV/c^2 applied in ROOT")

# J/psi veto: RM(pi+ pi-) < 3.09 or > 3.105 GeV/c^2 (veto psi(3686)->pi+pi-J/psi)
alg.note(:background_veto,
  "psi(3686)->pi+pi-J/psi veto: RM(pi+pi-) outside [3.09, 3.105] GeV/c^2; continuum background from 3.650/3.773 GeV data negligible")

# Lambda mass window [1.111, 1.120] and K_S0 mass window [0.489, 0.506] applied in ROOT
alg.note(:mass_windows,
  "Lambda mass window [1.111, 1.120]; K_S0 mass window [0.489, 0.506]; both 6-sigma windows applied in ROOT")

alg.with_decay_card(decay_card_signal).apply(event_selection)
alg.execute_on([psip_data, psip_incMC, exMC_signal])