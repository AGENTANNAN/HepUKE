# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")     # ψ(3686) real data (448.1×10^6 events)
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")  # corresponding inclusive MC sample

# Decay card for the full signal chain:
# ψ(3686) → π+π− J/ψ, J/ψ → φ e+e−, φ → K+K−
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.0000 pi+ pi- J/psi    PHSP;
    Enddecay

    Decay J/psi
    1.0000 phi e+ e-    PHOTOS VLL;
    Enddecay

    Decay phi
    1.0000 K+ K-    VSS;
    Enddecay

    End
DECAYCARD

# Exclusive MC sample for the full decay chain (100k events)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_pipi_jpsi_phi_ee"
  config.related_dataset = psip_data     # associated real dataset
  config.events          = 100000        # 100k signal MC events
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "PipiJpsiPhiEE"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})   # √s = 3.686 GeV
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection
    .select_track {                    # Charged track selection
        cos_theta 0.93                 # |cosθ| < 0.93
        Vz        10.0                 # |Vz| < 10 cm
        Vr        1.0                  # Vr < 1 cm
        nChrp     "==3"                # exactly 3 positive tracks
        nChrn     "==3"                # exactly 3 negative tracks
        nNet      "==0"                # net charge zero
    }
    .pid(method: :probability) {       # PID: probability method
        prob_cut 0.001                 # probability > 0.001
        # electrons/muons: high-momentum lepton identification (p>1.0 → lepton;
        # EMC energy > 0.6 → electron, else muon)
        identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                       treat_as_electron_if_energy_above: 0.6
        identify :kaon, against: [:pion, :proton]   # K+ and K− (both charges)
        nkp ">=1"                      # at least one K+
        nkm ">=1"                      # at least one K−
        nlp ">=1"                      # at least one e+
        nlm ">=1"                      # at least one e−
    }
    .remove([:kp <= :chrgp, :km <= :chrgn, :lp <= :chrgp, :lm <= :chrgn])  # remove identified K and e
    .assign({:chrgp => :pip, :chrgn => :pim})   # remaining tracks assigned as π+ / π−
    .remove(:pip) { condition "three_momentum_of(:pip) > 0.45" }   # soft pion: p < 0.45 GeV/c
    .remove(:pim) { condition "three_momentum_of(:pim) > 0.45" }   # soft pion: p < 0.45 GeV/c
    .recoil_mass_of(:pip, :pim).between(3.05, 3.15)   # M(π+π−)_rec selects the J/ψ
    .kinematic_fit([:pip, :pim, :kp, :km, :lp, :lm]) {  # 4C fit to π+π−K+K−e+e−
        nominal                     # mark as nominal fit
        constrain_four_momentum     # 4-momentum conservation against CMS energy
        chi2_cut 40                 # χ² < 40 (smallest-χ² assignment kept by default)
    }

# electron E/p cut cannot be expressed in the PID DSL block — capture it explicitly
my_algorithm
    .note(:electron_energy_momentum_cut,
          "electron candidates are additionally required to satisfy E/p > 0.8 (EMC energy over " \
          "MDC momentum); the PID DSL identifies high-momentum leptons via momentum and EMC energy " \
          "but does not express the E/p ratio, so this cut is applied in the generated BOSS algorithm")
    .with_decay_card(decay_card_signal)
    .apply(event_selection)

# Execute on real data, inclusive MC and the signal exclusive MC
root_files = my_algorithm.execute_on([psip_data, psip_incMC, exMC_signal])