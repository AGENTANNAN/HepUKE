# =============================================================================
# e+e- -> Lambda Lambdabar   (Lambda -> p pi-, Lambdabar -> anti-p pi+)
# Energy scan 3.51 - 4.60 GeV
# BOSS part: dataset preparation + event selection up to the 4C kinematic fit
# =============================================================================

### Dataset description ###
# Real data at each scan point of the XYZ energy scan (BOSS 703)
data_points = [
  DatasetManager.real_data.find("703_4180"),
  DatasetManager.real_data.find("703_4190"),
  DatasetManager.real_data.find("703_4200"),
  DatasetManager.real_data.find("703_4210"),
  DatasetManager.real_data.find("703_4220"),
  DatasetManager.real_data.find("703_4230"),
  DatasetManager.real_data.find("703_4260"),
  DatasetManager.real_data.find("703_4270"),
  DatasetManager.real_data.find("703_4280"),
  DatasetManager.real_data.find("703_4360"),
  DatasetManager.real_data.find("703_4420"),
  DatasetManager.real_data.find("703_4600")
]

# Corresponding inclusive MC samples at the same scan points
incMC_points = [
  DatasetManager.inclusive_mc.find("703_4180"),
  DatasetManager.inclusive_mc.find("703_4190"),
  DatasetManager.inclusive_mc.find("703_4200"),
  DatasetManager.inclusive_mc.find("703_4210"),
  DatasetManager.inclusive_mc.find("703_4220"),
  DatasetManager.inclusive_mc.find("703_4230"),
  DatasetManager.inclusive_mc.find("703_4260"),
  DatasetManager.inclusive_mc.find("703_4270"),
  DatasetManager.inclusive_mc.find("703_4280"),
  DatasetManager.inclusive_mc.find("703_4360"),
  DatasetManager.inclusive_mc.find("703_4420"),
  DatasetManager.inclusive_mc.find("703_4600")
]

# Decay card for the signal process (EvtGen format, EvtGen particle names).
# No intermediate resonance is produced directly, so the KKMC top mother
# psi(4260) is used as the incoming state.
decay_card_llbar = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Lambda0 anti-Lambda0            PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi-                          HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+                     HypWK;
    Enddecay

    End
DECAYCARD

# One 100k-event exclusive MC per scan point (same decay card, same cross section)
exMCs_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_llbar_scan"   # auto-suffixed per dataset, e.g. exmc_llbar_scan_703_4180
  config.events        = 100_000             # 100k events per scan point
  config.decay_card    = decay_card_llbar
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "LambdaLambdaBar"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 4.600]})  # nominal scan energy (GeV)

# BOSS-side procedures not expressible in the DSL are preserved as notes.
my_algorithm
  .note(:secondary_vertex_quality,
        "Lambda / Lambda_bar secondary-vertex candidates are required to have a
         vertex-fit chi2 < 500 and a positive decay length; these quality cuts are
         not expressible in the DSL and are applied on the virtual-particle
         candidates produced by secondary_vertex_fit.")
  .note(:lambda_mass_window,
        "Lambda mass window |M(p pi-) - m_Lambda| < 5 MeV/c^2 (and the charge
         conjugate for Lambda_bar) is required before the 4C fit; background is
         estimated from the Lambda mass sidebands (sideband subtraction) at the
         ROOT analysis level.")

# Build the event selection chain.
event_selection = Selection.new
event_selection
  .select_track {              # Charged track selection
      cos_theta  0.93          # |cos(theta)| < 0.93
      Vz         10.0          # |Vz| < 10 cm
      Vr         1.0           # Vr < 1 cm
      nChrp      "==2"         # Exactly two positive tracks
      nChrn      "==2"         # Exactly two negative tracks
      nNet       "==0"         # Net charge zero
  }
  .pid(method: :probability) { # Probability-based PID
      prob_cut   0.001         # PID probability > 0.001
      identify :proton, against: [:kaon, :pion]  # p+ and anti-p- vs K and pi
      nprp       ">=1"         # At least one proton
      nprm       ">=1"         # At least one anti-proton
  }
  .remove([:prp <= :chrgp])    # Remove identified protons from positive tracks
  .remove([:prm <= :chrgn])    # Remove identified anti-protons from negative tracks
  .assign({:chrgp => :pip, :chrgn => :pim})  # Remaining tracks treated as pi+ / pi-
  .secondary_vertex_fit([:prp, :pim]) {      # Lambda -> p pi-  secondary vertex
      build_virtual_particle(:Lambda).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:prm, :pip]) {      # Lambda_bar -> anti-p pi+  secondary vertex
      build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:Lambda, :Lambda_bar]) {   # 4C kinematic fit on the Lambda Lambda_bar system
      nominal                  # Nominal fit: its corrected four-momenta are used
      constrain_four_momentum  # Constrain total four-momentum to the CMS energy
      chi2_cut 200             # Loose chi2 < 200 (tight cut optimised in ROOT)
  }

# Attach the decay card and render the selection into the BOSS algorithm
my_algorithm.with_decay_card(decay_card_llbar).apply(event_selection)

# Execute on all scan points: real data + inclusive MC + signal exclusive MC
root_files = my_algorithm.execute_on(data_points + incMC_points + exMCs_signal)