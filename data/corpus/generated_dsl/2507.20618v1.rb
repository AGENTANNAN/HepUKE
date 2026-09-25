### Dataset description ###
psip_data = DatasetManager.real_data.find("709_3686")      # ψ(3686) real data at 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")  # Corresponding inclusive MC sample

# Decay card for the e+e- mode: ψ(3686) → π+π- J/ψ, J/ψ → e+e-
decay_card_signal_ee = <<~DECAYCARD
    Decay psi(2S)
    1  pi+  pi-  J/psi    PHSP;
    Enddecay

    Decay J/psi
    1  e+  e-    PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# Decay card for the μ+μ- mode: ψ(3686) → π+π- J/ψ, J/ψ → μ+μ-
decay_card_signal_mumu = <<~DECAYCARD
    Decay psi(2S)
    1  pi+  pi-  J/psi    PHSP;
    Enddecay

    Decay J/psi
    1  mu+  mu-    PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# Create 100k-event exclusive MC samples for each J/ψ decay mode
exMC_signal_ee = DatasetManager.create_exclusive_mc do |config|
    config.sample_name = "exmc_3686_pipijpsi_ee"
    config.related_dataset = psip_data
    config.events = 100000
    config.decay_card = decay_card_signal_ee
    config.cross_section = :default
end

exMC_signal_mumu = DatasetManager.create_exclusive_mc do |config|
    config.sample_name = "exmc_3686_pipijpsi_mumu"
    config.related_dataset = psip_data
    config.events = 100000
    config.decay_card = decay_card_signal_mumu
    config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "pipiJpsi"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})
            .set_alias({"std::vector<double>" => "Vdouble"})

# Common selection chain shared by both J/ψ decay modes
event_selection = Selection.new
event_selection.select_track {
                    cos_theta  0.93     # |cosθ| < 0.93
                    Vz         10.0     # |Vz| < 10 cm along the beam direction
                    Vr         1.0      # Vr < 1 cm in the transverse plane
                    nChrp      "==2"    # Exactly 2 positive tracks
                    nChrn      "==2"    # Exactly 2 negative tracks
                    nNet       "==0"    # Net charge zero
                }
                .pid(method: :probability) {
                    # High-momentum (p > 1.0 GeV/c) tracks treated as leptons;
                    # electron if EMC energy > 1.0 GeV, otherwise muon.
                    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                                   treat_as_electron_if_energy_above: 1.0
                    # Pions identified against kaons and protons
                    identify :pion, against: [:kaon, :proton]
                    npip  "==1"     # 1 π+
                    npim  "==1"     # 1 π-
                    nlp   "==1"     # 1 lepton l+  (e+ or μ+)
                    nlm   "==1"     # 1 lepton l-  (e- or μ-)
                }
                # 4C kinematic fit over π+π- l+l-
                .kinematic_fit([:pip, :pim, :lp, :lm]) do
                    nominal                              # Nominal fit — its corrected four-momenta are used
                    vertex_fit([0, 1, 2, 3])             # Vertex fit over all four tracks (pip, pim, lp, lm)
                    constrain_four_momentum              # 4C energy-momentum constraint
                    invariant_mass_of(:lp, :lm).within(3.087, 3.107)  # Constrain l+l- mass to J/ψ window
                    chi2_cut 200                         # χ² < 200
                end

# One common algorithm for both decay modes (identical final states and selection)
my_algorithm.with_decay_card(decay_card_signal_ee).apply(event_selection)

# Execute on real data, inclusive MC, and both exclusive MC samples
root_files = my_algorithm.execute_on([psip_data, psip_incMC, exMC_signal_ee, exMC_signal_mumu])