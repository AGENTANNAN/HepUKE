# =====================================================================
# Dataset preparation
# =====================================================================
jpsi_data  = DatasetManager.real_data.find("708_3097")     # J/psi real data at 3.097 GeV
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")  # corresponding inclusive MC

# Decay card for the antineutron tag: J/psi -> p+ pi- anti-n- (phase space).
# NOTE: the measured scattering final states K+ K- pi+ (pi0) originate from
# anti-n rescattering on 1H in the beam-pipe oil and are therefore not part
# of the generator-level card used for the tag exclusive MC.
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 p+ pi- anti-n- PHSP;
    Enddecay
    End
DECAYCARD

# 2M-event exclusive MC for J/psi -> p+ pi- anti-n-
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3097_ppimnbar"
  config.related_dataset = jpsi_data
  config.events          = 2_000_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

# =====================================================================
# Mode I : anti-n p -> K+ K- pi+   (reconstructed final state p pi- K+ K- pi+)
# =====================================================================
alg_name_I = "NbarPKKpi"
alg_I = Algorithm.new(alg_name_I)
alg_I.set_header(["#{alg_name_I}Alg/#{alg_name_I}.h"])
     .set_constant({"ECMS" => [:double, 3.097]})
     .set_alias({"std::vector<double>" => "Vdouble"})
     .note(:nbar_recoil_mass_constraint,
           "1C kinematic fit constraining the recoil mass against the p pi- system to the nominal "
           "antineutron mass (chi2 < 11); this 1C recoil-mass constraint has no dedicated DSL primitive.")
     .note(:nbar_tag_vertex_chi2,
           "p pi- secondary-vertex (anti-n tag) fit requires vertex chi2 < 5.")
     .note(:scattering_vertex_chi2,
           "K+ K- pi+ scattering-vertex fit requires vertex chi2 < 6.")
     .note(:impact_point_match,
           "scattering vertex must match the expected anti-n impact point on the beam pipe: "
           "|dr| < 0.3 cm and |dz| < 2.0 cm.")
     .note(:oil_layer_scattering,
           "anti-n scattering on 1H inside the oil layer is selected by the momentum balance "
           "P(p) = |P(K+ K- pi+) - P(anti-n)| < 0.04 GeV/c.")

sel_I = Selection.new
sel_I.select_track {
        cos_theta 0.93
        Vz        10.0
        Vr        1.0
        nChrp     ">=2"
        nChrn     ">=2"
      }
     .select_photon {
        tdc_emc_start     0
        tdc_emc_end       14
        angle_to_track    10.0
        energyThreshold_b 0.025
        energyThreshold_e 0.050
      }
     .pid(method: :probability) {
        prob_cut 0.001
        identify :proton, against: [:kaon, :pion]
        identify :kaon,   against: [:pion, :proton]
        identify :pion,   against: [:kaon, :proton]
        nprp ">=1"
        npim ">=1"
        nkp  ">=1"
        nkm  ">=1"
        npip ">=1"
      }
     # Tag the anti-n: secondary vertex fit of the p pi- system
     .secondary_vertex_fit([:prp, :pim]) {
        build_virtual_particle(:nbar_tag).by_minimizing_mass_difference
      }
     # 1C fit: recoil mass against p pi- constrained to the nominal anti-n mass
     .kinematic_fit([:prp, :pim]) {
        chi2_cut 11
      }
     # Scattering point: vertex fit of the K+ K- pi+ system
     .secondary_vertex_fit([:kp, :km, :pip]) {
        build_virtual_particle(:np_scatter).by_minimizing_verfit_chi2
      }
     # Final kinematic fit on p pi- K+ K- pi+
     .kinematic_fit([:prp, :pim, :kp, :km, :pip]) {
        nominal
        constrain_four_momentum
        chi2_cut 200
      }

alg_I.with_decay_card(decay_card_signal).apply(sel_I)

# =====================================================================
# Mode II : anti-n p -> K+ K- pi+ pi0  (final state p pi- K+ K- pi+ pi0)
# =====================================================================
alg_name_II = "NbarPKKpipi0"
alg_II = Algorithm.new(alg_name_II)
alg_II.set_header(["#{alg_name_II}Alg/#{alg_name_II}.h"])
      .set_constant({"ECMS" => [:double, 3.097]})
      .set_alias({"std::vector<double>" => "Vdouble"})
      .note(:nbar_recoil_mass_constraint,
            "1C kinematic fit constraining the recoil mass against the p pi- system to the nominal "
            "antineutron mass (chi2 < 11); this 1C recoil-mass constraint has no dedicated DSL primitive.")
      .note(:nbar_tag_vertex_chi2,
            "p pi- secondary-vertex (anti-n tag) fit requires vertex chi2 < 5.")
      .note(:scattering_vertex_chi2,
            "K+ K- pi+ scattering-vertex fit requires vertex chi2 < 6.")
      .note(:impact_point_match,
            "scattering vertex must match the expected anti-n impact point on the beam pipe: "
            "|dr| < 0.3 cm and |dz| < 2.0 cm.")
      .note(:oil_layer_scattering,
            "anti-n scattering on 1H inside the oil layer is selected by the momentum balance "
            "P(p) = |P(K+ K- pi+ pi0) - P(anti-n)| < 0.04 GeV/c.")

sel_II = Selection.new
sel_II.select_track {
        cos_theta 0.93
        Vz        10.0
        Vr        1.0
        nChrp     ">=2"
        nChrn     ">=2"
      }
      .select_photon {
        tdc_emc_start     0
        tdc_emc_end       14
        angle_to_track    10.0
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam              ">=2"
      }
      .pid(method: :probability) {
        prob_cut 0.001
        identify :proton, against: [:kaon, :pion]
        identify :kaon,   against: [:pion, :proton]
        identify :pion,   against: [:kaon, :proton]
        nprp ">=1"
        npim ">=1"
        nkp  ">=1"
        nkm  ">=1"
        npip ">=1"
      }
      # pi0 reconstruction from photon pairs (1C Kalman fit)
      .kalman_kinematic_fit([:gamma, :gamma]) {
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 50
        npi0 ">=1"
      }
      # Tag the anti-n: secondary vertex fit of the p pi- system
      .secondary_vertex_fit([:prp, :pim]) {
        build_virtual_particle(:nbar_tag).by_minimizing_mass_difference
      }
      # 1C fit: recoil mass against p pi- constrained to the nominal anti-n mass
      .kinematic_fit([:prp, :pim]) {
        chi2_cut 11
      }
      # Scattering point: vertex fit of the K+ K- pi+ system
      .secondary_vertex_fit([:kp, :km, :pip]) {
        build_virtual_particle(:np_scatter).by_minimizing_verfit_chi2
      }
      # Final kinematic fit on p pi- K+ K- pi+ pi0
      .kinematic_fit([:prp, :pim, :kp, :km, :pip, :pi0]) {
        nominal
        constrain_four_momentum
        chi2_cut 200
      }

alg_II.with_decay_card(decay_card_signal).apply(sel_II)

# =====================================================================
# Execution on data / inclusive MC / signal exclusive MC
# =====================================================================
root_files_I  = alg_I.execute_on([jpsi_data, jpsi_incMC, exMC_signal])
root_files_II = alg_II.execute_on([jpsi_data, jpsi_incMC, exMC_signal])