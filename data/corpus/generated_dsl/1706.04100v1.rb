### Dataset description ###
# Real data and inclusive MC at the two ψ(4260) energy points (4.226 and 4.258 GeV)
data_4230  = DatasetManager.real_data.find("703_4230")      # Ecms ≈ 4226 MeV (4.226 GeV)
data_4260  = DatasetManager.real_data.find("703_4260")      # Ecms ≈ 4258 MeV (4.258 GeV)
incMC_4230 = DatasetManager.inclusive_mc.find("703_4230")   # corresponding inclusive MC
incMC_4260 = DatasetManager.inclusive_mc.find("703_4260")

# Decay card: ψ(4260) → π⁺π⁻ J/ψ, J/ψ → e⁺e⁻ in phase space
decay_card_ee = <<~DECAYCARD
  Decay psi(4260)
  1.000  pi+  pi-  J/psi    PHSP;
  Enddecay

  Decay J/psi
  1.000  e+  e-    PHSP;
  Enddecay

  End
DECAYCARD

# Decay card: ψ(4260) → π⁺π⁻ J/ψ, J/ψ → μ⁺μ⁻ in phase space
decay_card_mumu = <<~DECAYCARD
  Decay psi(4260)
  1.000  pi+  pi-  J/psi    PHSP;
  Enddecay

  Decay J/psi
  1.000  mu+  mu-    PHSP;
  Enddecay

  End
DECAYCARD

# 100k-event exclusive MC for each J/ψ decay mode at both energy points
exMC_ee = DatasetManager.create_exclusive_mc_for([data_4230, data_4260]) do |config|
  config.sample_name   = "exmc_pipijpsi_ee"
  config.events        = 100_000
  config.decay_card    = decay_card_ee
  config.cross_section = :default
end

exMC_mumu = DatasetManager.create_exclusive_mc_for([data_4230, data_4260]) do |config|
  config.sample_name   = "exmc_pipijpsi_mumu"
  config.events        = 100_000
  config.decay_card    = decay_card_mumu
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "pipiJpsi"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 4.226]})
            .note(:energy_point, "selection shared by the 4.226 and 4.258 GeV data points; a single ECMS constant cannot cover both, so the per-dataset beam energy is used for the kinematic fit at execute_on")

# Both J/ψ decay modes (e⁺e⁻ and μ⁺μ⁻) share this single selection chain
event_selection = Selection.new
  .select_track {
    cos_theta 0.93       # |cosθ| < 0.93
    Vz        10.0       # |Vz| < 10 cm
    Vr        1.0        # Vr < 1 cm
    nChrp     "==2"      # exactly two positive tracks
    nChrn     "==2"      # exactly two negative tracks
    nNet      "==0"      # net charge zero
  }
  .pid(method: :probability) {
    # High-momentum tracks (p > 1.0) treated as leptons; electron if EMC energy > 0.6, else muon
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6
    identify :pion, against: [:kaon]   # π⁺ / π⁻ from ψ(4260) → π⁺π⁻ J/ψ
    npip "==1"
    npim "==1"
    nlp  "==1"                          # J/ψ → l⁺ l⁻
    nlm  "==1"
  }
  # 5C kinematic fit: 4C (e⁺e⁻ → π⁺π⁻ J/ψ) + J/ψ mass constraint on the lepton pair
  .kinematic_fit([:pip, :pim, :lp, :lm]) {
    nominal
    constrain_four_momentum
    invariant_mass_of(:lp, :lm).constrain_to_nominal_mass_of(:jpsi)
    invariant_mass_of(:lp, :lm).within(0.1, 3.2)   # J/ψ mass window on M(l⁺l⁻)
    chi2_cut 200
  }

my_algorithm.with_decay_card(decay_card_ee).apply(event_selection)

root_files = my_algorithm.execute_on(
  [data_4230, data_4260, incMC_4230, incMC_4260] + exMC_ee + exMC_mumu
)