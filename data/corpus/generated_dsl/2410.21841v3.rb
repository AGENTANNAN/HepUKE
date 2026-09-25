# =============================================================================
# J/psi (sqrt(s) = 3.097 GeV) search for Lambda-Lambdabar oscillation
#   right-sign : J/psi -> Lambda anti-Lambda
#   wrong-sign : J/psi -> Lambda Lambda      (oscillation mode)
#   peaking bg : J/psi -> Lambda anti-Sigma0 + c.c.
# =============================================================================

### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")       # J/psi real data at 3.097 GeV
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")    # Corresponding inclusive MC

# --- Decay card: right-sign signal J/psi -> Lambda anti-Lambda ---
decay_card_rs = <<~DECAYCARD
    Decay J/psi
    1.0000 Lambda0 anti-Lambda0 PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+ HypWK;
    Enddecay

    End
DECAYCARD

# --- Decay card: wrong-sign (oscillation) signal J/psi -> Lambda Lambda ---
decay_card_ws = <<~DECAYCARD
    Decay J/psi
    1.0000 Lambda0 Lambda0 PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- HypWK;
    Enddecay

    End
DECAYCARD

# --- Decay card: peaking background J/psi -> Lambda anti-Sigma0 + c.c. ---
decay_card_bg = <<~DECAYCARD
    Decay J/psi
    0.5000 Lambda0 anti-Sigma0 PHSP;
    0.5000 anti-Lambda0 Sigma0 PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+ HypWK;
    Enddecay

    Decay anti-Sigma0
    1.0000 gamma anti-Lambda0 PHSP;
    Enddecay

    Decay Sigma0
    1.0000 gamma Lambda0 PHSP;
    Enddecay

    End
DECAYCARD

# --- Exclusive MC samples ---
exMC_rs = DatasetManager.create_exclusive_mc do |config|   # 10M right-sign J/psi -> Lambda anti-Lambda
  config.sample_name     = "exmc_jpsi_LLbar"
  config.related_dataset = jpsi_data
  config.events          = 10_000_000
  config.decay_card      = decay_card_rs
  config.cross_section   = :default
end

exMC_ws = DatasetManager.create_exclusive_mc do |config|   # 1M wrong-sign (oscillation) J/psi -> Lambda Lambda
  config.sample_name     = "exmc_jpsi_LL_osc"
  config.related_dataset = jpsi_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_ws
  config.cross_section   = :default
end

exMC_bg = DatasetManager.create_exclusive_mc do |config|   # 0.5M peaking background J/psi -> Lambda anti-Sigma0 + c.c.
  config.sample_name     = "exmc_jpsi_LSigmabar0"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_bg
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "LambdaLambdaBarOsc"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})   # sqrt(s) = 3.097 GeV
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {          # Charged track selection
                  cos_theta 0.93        # |cos(theta)| < 0.93
                  Vz        10.0        # |Vz| < 10 cm
                  Vr        1.0         # Vr < 1 cm
                  nTot      ">=4"       # At least four charged tracks
                }
               .pid(method: :probability) {   # Probability-based PID
                  prob_cut 0.001               # PID probability > 0.001
                  identify :proton, against: [:kaon, :pion]   # p+ and anti-p- vs K and pi
                  nprp ">=1"                   # At least one proton
                  nprm ">=1"                   # At least one anti-proton
                }
               .remove([:prp <= :chrgp])       # Remove identified protons from positive list
               .remove([:prm <= :chrgn])       # Remove identified anti-protons from negative list
               .assign({:chrgp => :pip, :chrgn => :pim})   # Remaining tracks taken as pi+ and pi-
               .secondary_vertex_fit([:prp, :pim]) {   # Reconstruct Lambda from p pi-
                    build_virtual_particle(:Lambda).by_minimizing_mass_difference
                    remove_used_particle_from_candidate_list
                }
               .secondary_vertex_fit([:prm, :pip]) {   # Reconstruct anti-Lambda from anti-p pi+
                    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
                    remove_used_particle_from_candidate_list
                }
               .invariant_mass_of(:Lambda).between(1.1103, 1.1211)      # +-3 sigma window (sigma = 1.8 MeV)
               .invariant_mass_of(:Lambda_bar).between(1.1103, 1.1211)  # +-3 sigma window (sigma = 1.8 MeV)
               .kinematic_fit([:Lambda, :Lambda_bar]) {   # Nominal 4C kinematic fit to Lambda anti-Lambda
                    nominal                  # best-chi2 combination is chosen automatically
                    constrain_four_momentum  # 4C energy-momentum conservation
                    chi2_cut 50              # chi2 < 50
                }

my_algorithm.with_decay_card(decay_card_rs).apply(event_selection)

root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_rs, exMC_ws, exMC_bg])