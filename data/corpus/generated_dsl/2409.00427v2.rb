# =====================================================================
# BOSS DSL specification
# Born cross section of e+e- -> Xi0 anti-Xi0 (single-baryon-tag method,
# charge conjugate included).
# Decay chain: Xi0 -> pi0 Lambda0, Lambda0 -> p+ pi-, pi0 -> gamma gamma.
# Scan: 45 CM-energy points, 3.51-4.95 GeV  (see assets/BES3_dataset.md).
# =====================================================================

### ------------------------------------------------------------------ ###
### Dataset preparation                                                ###
### ------------------------------------------------------------------ ###

# The 45 CM-energy scan points, sample-name convention "<boss>_<ecms_MeV>"
scan_point_names = [
  # BOSS 703
  "703_3810", "703_3872", "703_3900", "703_4009", "703_4090",
  "703_4180", "703_4190", "703_4200", "703_4210", "703_4220",
  "703_4230", "703_4237", "703_4245", "703_4246", "703_4260",
  "703_4270", "703_4280", "703_4310", "703_4360", "703_4390",
  "703_4420", "703_4470", "703_4530", "703_4575", "703_4600",
  # BOSS 705
  "705_4130", "705_4160", "705_4290", "705_4315", "705_4340",
  "705_4380", "705_4400", "705_4440",
  # BOSS 706
  "706_4610", "706_4620", "706_4640", "706_4660", "706_4680", "706_4700",
  # BOSS 707
  "707_4740", "707_4750", "707_4780", "707_4840", "707_4914", "707_4946"
]

scan_data  = scan_point_names.map { |n| DatasetManager.real_data.find(n) }     # real data
scan_incMC = scan_point_names.map { |n| DatasetManager.inclusive_mc.find(n) } # matching inclusive MC

# ---- Decay card: continuum production, psi(4260) top mother, flat phase space ----
decay_card_signal = <<~DECAYCARD
  Decay psi(4260)
  1.000 Xi0 anti-Xi0 PHSP;
  Enddecay

  Decay Xi0
  1.000 pi0 Lambda0 PHSP;
  Enddecay

  Decay anti-Xi0
  1.000 pi0 anti-Lambda0 PHSP;
  Enddecay

  Decay Lambda0
  1.000 p+ pi- HypWK;
  Enddecay

  Decay anti-Lambda0
  1.000 anti-p- pi+ HypWK;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# ---- Exclusive signal MC: 100k events at every scan point ----
exMC_signal = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name   = "exmc_xi0xibar0_scan"   # auto-suffixed per energy point
  config.events        = 100_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end
exMC_signal.each { |m| m.save_to_config(format: :yaml, file_path: 'temp_for_test') }

### ------------------------------------------------------------------ ###
### Event selection (BOSS)                                            ###
### ------------------------------------------------------------------ ###

alg_name = "Xi0Xi0barST"
xi0_alg = Algorithm.new(alg_name)
xi0_alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
       .set_constant({"ECMS" => [:double, 4.260]})  # nominal value; varies over the scan (see note)
       .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
  .select_track {                  # charged-track selection
    cos_theta 0.93                 # |cos(theta)| < 0.93
    Vz        100.0                # |Vz| < 100 cm
    Vr        10.0                 # Vr < 10 mm
    nChrp     ">=2"                # at least two positive tracks
    nChrn     ">=1"                # at least one negative track
    nNet      "==0"                # net charge zero
  }
  .select_photon {                 # photon selection
    tdc_emc_start     0            # EMC timing 0-14
    tdc_emc_end       14
    angle_to_track    10.0         # at least 10 deg from any charged track
    energyThreshold_b 0.025        # 25 MeV (barrel)
    energyThreshold_e 0.050        # 50 MeV (endcap)
    nGam              ">=2"        # at least two photons
  }
  .pid(method: :probability) {     # PID: protons against kaons and pions
    prob_cut 0.001                 # probability cut
    identify :proton, against: [:kaon, :pion]   # p+ and anti-p- at once
    nprp ">=1"                     # at least one proton
    nprm ">=1"                     # at least one antiproton
  }
  .remove([:prp <= :chrgp, :prm <= :chrgn])   # remove proton tracks from charged lists
  .assign({:chrgp => :pip, :chrgn => :pim})   # remaining + -> pi+, - -> pi-
  .secondary_vertex_fit([:prp, :pim]) {       # Lambda0 -> p+ pi-
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:prm, :pip]) {       # anti-Lambda0 -> anti-p- pi+
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # pi0 -> gamma gamma (1C mass constraint)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 20
    npi0 ">=1"                                # at least one pi0
  }
  # No kinematic fit is applied to the Xi0: the tag-side Xi0 (recID 1) is built
  # from pi0 + Lambda0 and the untagged anti-Xi0 is inferred from the recoil.
  .partial_rec([1]) {
    best_combination_by_mass :Xi0, 1.31486    # minimise |M(pi0 Lambda) - m_Xi0|
  }

# Inexpressible / downstream (ROOT) procedures preserved as notes
xi0_alg
  .note(:root_level_xi0_selection,
    "No kinematic fit is applied to the Xi0 because it is not directly measured. " \
    "The Xi0 candidate is chosen by minimising |M(pi0 Lambda) - m_Xi0|; the signal " \
    "window |M(pi0 Lambda) - m_Xi0| < 10 MeV/c2, the recoil-mass window " \
    "|M_recoil - m_Xi0| < 60 MeV/c2 and the Lambda0 window |M(p pi-) - m_Lambda| < 5 MeV/c2 " \
    "are all applied in the downstream ROOT stage.")
  .note(:background_estimate,
    "Four-region two-dimensional sideband background estimate applied in the " \
    "downstream ROOT stage.")
  .note(:charge_conjugation,
    "Charge conjugate included: the same event topology with the anti-Xi0 tag " \
    "(anti-Xi0 -> pi0 anti-Lambda0) is covered by the identical selection.")
  .note(:ecms_scan,
    "ECMS is not one constant: the 45 scan points span 3.51-4.95 GeV. The value set " \
    "here is nominal and must be overridden per energy point at execution.")
  .with_decay_card(decay_card_signal)
  .apply(event_selection)

# Execute on real data, matching inclusive MC and the scan-point signal MC
root_files = xi0_alg.execute_on(scan_data + scan_incMC + exMC_signal)