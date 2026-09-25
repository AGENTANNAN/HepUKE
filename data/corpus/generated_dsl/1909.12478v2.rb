# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
# Seven energy points of the BESIII e+e- scan (BOSS 7.0.3):
# 4358.3, 4387.4, 4415.6, 4467.1, 4527.1, 4574.5, 4599.5 MeV
energy_points = %w[703_4360 703_4390 703_4420 703_4470 703_4530 703_4575 703_4600]
data_samples  = energy_points.map { |name| DatasetManager.real_data.find(name) }
incMC_samples = energy_points.map { |name| DatasetManager.inclusive_mc.find(name) }

# Decay card for the signal process e+e- -> D+ D- pi+ pi- (EvtGen format).
# The same final state is populated both by e+e- -> D1(2420)+ D- + c.c.
# (D1(2420)+ -> D+ pi+ pi-) and by e+e- -> psi(3770) pi+ pi- (psi(3770) -> D+ D-);
# only the D+ -> K- pi+ pi+ decay is reconstructed, the D- is taken as recoil.
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0000 D+ D- pi+ pi- PHSP;
    Enddecay

    Decay D+
    1.0000 K- pi+ pi+ PHSP;
    Enddecay

    Decay D-
    1.0000 K+ pi- pi- PHSP;
    Enddecay

    End
DECAYCARD

# 100k-event exclusive (signal) MC for the D+ D- pi+ pi- final state at each energy point
exMCs_signal = DatasetManager.create_exclusive_mc_for(data_samples) do |config|
    config.sample_name   = "exmc_DDbarPiPi"   # auto-suffixed per energy point
    config.events        = 100_000
    config.decay_card    = decay_card_signal
    config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "D1D"
alg = Algorithm.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 4.4156]})  # representative scan energy; set per point when running
   .note(:vertex_fit_chi2, "the D+ -> K- pi+ pi+ secondary vertex fit is required to satisfy chi2 < 100; secondary_vertex_fit exposes no chi2 cut, so the best combination is chosen by minimising the mass difference and the chi2 < 100 requirement is imposed in the generated BOSS code")
   .note(:missing_d_minus, "the D- is not reconstructed: it is inferred from the recoil mass against the D+ pi+ pi- system, and the combinatorial background is removed by sideband subtraction. The kinematic fit constrains only the visible D+ pi+ pi- system")
   .note(:energy_scan, "seven-point scan (4.358-4.600 GeV); the ECMS constant must be set to the corresponding value for each energy point")

event_selection = Selection.new
    .select_track {                        # charged track selection
        cos_theta 0.93                     # |cos(theta)| < 0.93
        Vz        100.0                    # |Vz| < 100 cm
        Vr        10.0                     # Vr < 10 cm
        nChrp     ">=3"                    # at least 3 positive tracks
        nChrn     ">=3"                    # at least 3 negative tracks
    }
    .pid(method: :probability) {           # PID with the probability method
        prob_cut 0.001                     # PID probability > 0.001
        identify :kaon, against: [:pion, :proton]  # K+ / K- separated from pi and p
        nkm ">=1"                          # at least one K-
    }
    .remove([:kp <= :chrgp, :km <= :chrgn])     # do not re-count identified kaons as pions
    .assign({:chrgp => :pip, :chrgn => :pim})   # remaining tracks assigned as pi+ / pi-
    .secondary_vertex_fit([:km, :pip, :pip]) {  # vertex fit for D+ -> K- pi+ pi+
        build_virtual_particle(:Dplus).by_minimizing_mass_difference
    }
    .kinematic_fit([:km, :pip, :pip, :pip, :pim]) {  # D+ mass constraint + 4-momentum conservation
        nominal
        invariant_mass_of(:km, :pip, :pip).constrain_to_nominal_mass_of(:Dplus)
        constrain_four_momentum
        chi2_cut 20
    }

# Generate the complete algorithm for the process defined by the decay card
alg.with_decay_card(decay_card_signal).apply(event_selection)

# Execute on real data, inclusive MC and the signal MC at all seven energy points
root_files = alg.execute_on(data_samples + incMC_samples + exMCs_signal)