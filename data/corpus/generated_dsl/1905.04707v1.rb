### Dataset description ###
data_4600  = DatasetManager.real_data.find("703_4600")        # 4.600 GeV data (~567 pb^-1)
incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")     # Matching inclusive MC at 4.600 GeV

# Decay card for e+e- -> Lambda_c+ anti-Lambda_c- (partial reconstruction of one Lambda_c+)
# Top mother psi(4260) per BESIII KKMC convention. Lambda_c+ -> p K_S0, K_S0 -> pi+ pi-.
# The anti-Lambda_c- side is left undecayed so it stays as the single missed recoil particle.
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0 Lambda_c+ anti-Lambda_c- PHSP;
    Enddecay

    Decay Lambda_c+
    1.0 p+ K_S0 PHSP;
    Enddecay

    Decay K_S0
    1.0 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# 500k exclusive MC events for Lambda_c+ -> p K_S0 (K_S0 -> pi+ pi-); anti-Lambda_c- not reconstructed
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_4600_LcToPKs"
  config.related_dataset = data_4600
  config.events          = 500000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "LcToPKsPartRec"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 4.600]})          # CMS energy 4.600 GeV
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {                 # Charged track selection
                  cos_theta 0.93               # |cos(theta)| < 0.93
                  Vz        10.0               # |Vz| < 10 cm
                  Vr        1.0                # Vr < 1 cm
                  nChrp     ">=2"              # at least 2 positive tracks
                  nChrn     ">=1"              # at least 1 negative track
                  nNet      "==1"              # net charge +1
                }
               # No photon requirement is imposed.
               .pid(method: :probability) {    # Particle identification (probability method)
                  prob_cut 0.001               # PID probability > 0.001
                  identify :proton, against: [:pion, :kaon]   # p+ / anti-p- vs pi, K
                  nprp ">=1"                   # at least one proton
               }
               .remove([:prp <= :chrgp])       # remove the identified proton from the charged list
               .assign({:chrgp => :pip, :chrgn => :pim})   # remaining positive/negative tracks -> pi+ / pi-
               .secondary_vertex_fit([:pip, :pim]) {       # form K_S0 -> pi+ pi-
                  build_virtual_particle(:K_S0).by_minimizing_mass_difference  # best candidate by mass difference
                  remove_used_particle_from_candidate_list
               }
               .partial_miss([2]) {            # recID 2 = anti-Lambda_c- (missed); no kinematic fit is applied
                  best_combination_by_mass :Lambda_c, 2.28646   # best p + K_S0 combination closest to nominal Lambda_c+ mass
                  require_recoil_mass 2.278, 2.294               # M_BC window on recoil mass against the missed anti-Lambda_c-
               }

my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)
root_files = my_algorithm.execute_on([data_4600, incMC_4600, exMC_signal])