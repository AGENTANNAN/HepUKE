# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset preparation ###
# J/psi line-shape scan: 15 centre-of-mass points, 3049.642 - 3119.878 MeV
scan_energies = [3049.642, 3054.659, 3059.676, 3064.693, 3069.709,
                 3074.726, 3079.743, 3084.760, 3089.777, 3094.794,
                 3099.810, 3104.827, 3109.844, 3114.861, 3119.878]

# Real data and inclusive MC (background estimation) at every scan point
scan_data  = scan_energies.map { |e| DatasetManager.real_data.find("708_#{e.round}") }
scan_incMC = scan_energies.map { |e| DatasetManager.inclusive_mc.find("708_#{e.round}") }

# ConExc decay cards for the two leptonic final states (continuum / ISR Born cross
# section). The DSL auto-detects the literal "ConExc" token, switches to the no-KKMC
# template and injects "Particle vpho <ECMS> 0.0" per energy point, so no explicit
# Particle vpho line is written for this multi-energy scan.
decay_card_bhabha = <<~DECAYCARD
    Decay vpho
    1.0 ConExc 1;
    Enddecay
    End
DECAYCARD

decay_card_dimuon = <<~DECAYCARD
    Decay vpho
    1.0 ConExc 2;
    Enddecay
    End
DECAYCARD

# 500k-event exclusive MC for each of the two leptonic final states, generated with
# the ConExc card and one sample per scan point.
exMC_bhabha = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name   = "jpsi_scan_bhabha"
  config.events        = 500_000
  config.decay_card    = decay_card_bhabha
  config.cross_section = :default
end

exMC_dimuon = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name   = "jpsi_scan_dimuon"
  config.events        = 500_000
  config.decay_card    = decay_card_dimuon
  config.cross_section = :default
end

### Bhabha channel: e+e- -> e+e- ###
bhabha_name = "JpsiBhabha"
bhabha_alg  = Algorithm.new(bhabha_name)
bhabha_alg.set_header(["#{bhabha_name}Alg/#{bhabha_name}.h"])
          .set_constant({"ECMS" => [:double, 3.085]})   # scan midpoint; per-point beam energy comes from the data
          .set_alias({"std::vector<double>" => "Vdouble"})
          .note(:lepton_selection,
                "Per-track Bhabha selection P > 0.7 E_beam and E > 0.6 P (E/p electron identification) has no dedicated DSL method; applied in the generated BOSS algorithm with E_beam = ECMS/2.")
          .note(:isr_vp_factor,
                "The ISR / vacuum-polarisation correction factor f is read from the ConExc generator log at each scan point and enters sigma = (N_sig - N_bkg)/(L * eps_trg * eps_sel) * f.")

bhabha_selection = Selection.new
bhabha_selection.select_track {
        cos_theta 0.8    # |cos(theta)| < 0.8 for both tracks
        Vz 10.0          # |Vz| < 10 cm
        Vr 1.0           # Vr < 1 cm
        nChrp "==1"      # exactly one positive track
        nChrn "==1"      # exactly one negative track
        nNet  "==0"      # net charge zero -> two oppositely charged tracks
      }
      .assign({:chrgp => :ep, :chrgn => :em})   # the E/p cut identifies the two electrons
      .kinematic_fit([:ep, :em]) {
        nominal
        constrain_four_momentum
        chi2_cut 200
      }

bhabha_alg.with_decay_card(decay_card_bhabha).apply(bhabha_selection)

### Dimuon channel: e+e- -> mu+mu- ###
dimuon_name = "JpsiDimuon"
dimuon_alg  = Algorithm.new(dimuon_name)
dimuon_alg.set_header(["#{dimuon_name}Alg/#{dimuon_name}.h"])
          .set_constant({"ECMS" => [:double, 3.085]})
          .set_alias({"std::vector<double>" => "Vdouble"})
          .note(:lepton_selection,
                "Per-track dimuon selection P > 0.8 E_beam, 25 MeV < E < 0.25 P and |Delta t_TOF| < 1.5 ns (cosmic suppression) has no dedicated DSL method; applied in the generated BOSS algorithm with E_beam = ECMS/2.")
          .note(:isr_vp_factor,
                "The ISR / vacuum-polarisation correction factor f is read from the ConExc generator log at each scan point and enters sigma = (N_sig - N_bkg)/(L * eps_trg * eps_sel) * f.")

dimuon_selection = Selection.new
dimuon_selection.select_track {
        cos_theta 0.8    # |cos(theta)| < 0.8 for both tracks
        Vz 10.0          # |Vz| < 10 cm
        Vr 1.0           # Vr < 1 cm
        nChrp "==1"      # exactly one positive track
        nChrn "==1"      # exactly one negative track
        nNet  "==0"      # net charge zero -> two oppositely charged tracks
      }
      .select_photon {
        tdc_emc_start 0
        tdc_emc_end 14
        angle_to_track 10.0
        energyThreshold_b 0.025
        energyThreshold_e 0.025
        nGam "==0"       # no neutral shower above 25 MeV
      }
      .assign({:chrgp => :mup, :chrgn => :mum})   # the energy / E-p cuts identify the two muons
      .kinematic_fit([:mup, :mum]) {
        nominal
        constrain_four_momentum
        chi2_cut 200
      }

dimuon_alg.with_decay_card(decay_card_dimuon).apply(dimuon_selection)

### Execution ###
# Each channel is run on all 15 data points + inclusive MC + its own exclusive MC
bhabha_files = bhabha_alg.execute_on(scan_data + scan_incMC + exMC_bhabha)
dimuon_files = dimuon_alg.execute_on(scan_data + scan_incMC + exMC_dimuon)