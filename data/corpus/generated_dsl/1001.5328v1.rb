### Dataset description ###
psip_data = DatasetManager.real_data.find("709_3686")      # ψ(2S) real data at 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")  # corresponding inclusive MC sample

# Decay card for ψ(2S) → π+π− J/ψ, J/ψ → γ p p̄ (phase space)
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.0000 pi+ pi- J/psi    PHSP;
    Enddecay

    Decay J/psi
    1.0000 gamma p+ anti-p-  PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC for the signal process (200k events)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3686_pipijpsi_gammappbar"
  config.related_dataset = psip_data
  config.events          = 200000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "PsipPipiJpsiGammaPPbar"
my_Algorithm = Algorithm.new(alg_name)
my_Algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {          # Charged track selection
                  cos_theta 0.93        # |cosθ| < 0.93
                  Vz        100.0       # |Vz| < 100 cm
                  Vr        10.0        # Vr < 10 mm
                  nChrp     ">=2"       # at least two positive tracks
                  nChrn     ">=2"       # at least two negative tracks
                  nNet      "==0"       # net charge zero
                }
               .select_photon {         # Photon selection
                  tdc_emc_start     0   # EMC TDC start
                  tdc_emc_end       14  # EMC TDC end
                  energyThreshold_b 0.025  # > 25 MeV in the barrel
                  energyThreshold_e 0.050  # > 50 MeV in the endcap
                  angle_to_track    10.0   # isolated from all charged tracks by > 10°
                  nGam              ">=1"  # at least one photon
                }
               .pid(method: :probability) {   # Particle identification
                  prob_cut 0.001                  # PID probability > 0.001
                  identify :proton, against: [:kaon, :pion]  # identify p+ and p̄ vs K, π
                  nprp ">=1"                      # at least one proton
                  nprm ">=1"                      # at least one antiproton
                }
               .remove([:prp <= :chrgp, :prm <= :chrgn])  # remove identified p / p̄ from the charged lists
               .assign({:chrgp => :pip, :chrgn => :pim})   # remaining tracks taken as π+ / π−
               .select_isolated_photon {       # Stricter isolation from the antiproton
                  angle_to_prm_track 30.0      # photon isolated from p̄ by > 30°
                  nGam ">=1"                   # at least one surviving photon
                }
               .kinematic_fit([:gamma, :pip, :pim, :prp, :prm]) {  # 4C fit to γ π+π− p p̄
                  nominal                      # nominal fit (corrected four-momenta saved)
                  constrain_four_momentum      # four-momentum constraint to the CMS
                  chi2_cut 100                 # χ² < 100
                }

# J/ψ-tag recoil-mass window on the π+π− system (|M_recoil(π+π−) − m_J/ψ| < ...):
# not expressible with current DSL constructs (recoil_mass_of is not implemented) and
# the bound was truncated in the source description — captured as a note.
my_Algorithm.note(:mrec_pipi_jpsi_window,
  "J/psi-tag recoil-mass window applied to the π+π− system: require the recoil mass
   against π+π− to be consistent with the J/psi nominal mass (|M_recoil(π+π−) − m_J/ψ|
   < bound). The bound value was truncated in the source description. No DSL primitive
   exists for the recoil mass of a subsystem (recoil_mass_of not implemented).")

my_Algorithm.with_decay_card(decay_card_signal).apply(event_selection)

root_files = my_Algorithm.execute_on([psip_data, psip_incMC, exMC_signal])