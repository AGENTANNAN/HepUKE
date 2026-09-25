### Dataset description ###
psip_data = DatasetManager.real_data.find("709_3686")        # ψ(3686) real data (BOSS + CMS energy convention)
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")    # Corresponding inclusive MC sample for ψ(3686)

# Decay card for the signal process ψ(2S) → π+π− J/ψ, J/ψ → γη, η → γγ
# (final state π+π−γγγ: two photons from η, one photon from J/ψ)
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.0000 pi+ pi- J/psi    PHSP;
    Enddecay

    Decay J/psi
    1.0000 gamma eta    PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma    PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for the J/ψ → γπ0π0 peaking background (estimated from dedicated MC)
decay_card_bkg = <<~DECAYCARD
    Decay psi(2S)
    1.0000 pi+ pi- J/psi    PHSP;
    Enddecay

    Decay J/psi
    1.0000 gamma pi0 pi0    PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma    PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC for the signal process (500k events)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
    config.sample_name = "exmc_psip_pipijpsi_gammaeta"
    config.related_dataset = psip_data      # Associated real dataset
    config.events = 500000                  # 500k signal events
    config.decay_card = decay_card_signal
    config.cross_section = :default
end

# Exclusive MC for the peaking background J/ψ → γπ0π0
exMC_bkg = DatasetManager.create_exclusive_mc do |config|
    config.sample_name = "exmc_psip_pipijpsi_gammapi0pi0"
    config.related_dataset = psip_data
    config.events = 100000
    config.decay_card = decay_card_bkg
    config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "PsipGammaEta"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})   # ECMS = 3.686 GeV

event_selection = Selection.new
event_selection.select_track {                # Charged track selection
                  cos_theta 0.93              # |cosθ| < 0.93
                  Vz        10.0              # |Vz| < 10 cm
                  Vr        1.0               # Vr < 1 cm
                  nChrp     "==1"             # Exactly one positive track
                  nChrn     "==1"             # Exactly one negative track
                  nNet      "==0"             # Net charge zero
                }
               .assign({:chrgp => :pip, :chrgn => :pim})   # Assume π+ and π−, no further PID
               .select_photon {               # Photon selection
                  tdc_emc_start     0         # EMC TDC start
                  tdc_emc_end       14        # EMC TDC end
                  angle_to_track    15.0      # > 15° from any charged track
                  energyThreshold_b 0.025     # > 25 MeV in EMC barrel
                  energyThreshold_e 0.050     # > 50 MeV in EMC endcap
                  nGam              ">=3"     # At least three photons
                }
               .kinematic_fit([:pip, :pim, :gamma, :gamma, :gamma]) {  # 4C fit to π+π−γγγ
                  nominal                         # Nominal kinematic fit
                  constrain_four_momentum         # Constrain total four-momentum to CMS energy
                  chi2_cut 50                     # χ² < 50
                }

# Attach signal decay card, generate algorithm, execute on all datasets
my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)
root_files = my_algorithm.execute_on([psip_data, psip_incMC, exMC_signal, exMC_bkg])