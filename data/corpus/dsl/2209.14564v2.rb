# Paper: 2209.14564v2
# Title: Observation of psi(3686) -> Sigma- Sigma+bar and measurement of its angular distribution
# Energy: 3.686 GeV (psi(3686), single energy point)
# Final state: psi(3686) -> Sigma- Sigma+bar; Sigma- -> n pi-; Sigma+bar -> nbar pi+
# Visible: pi+ pi- + anti-neutron (reconstructed via EMC showers)
# Neutron treated as missing particle (not reconstructed)

### Dataset preparation ###
data_709_3686 = DatasetManager.load_real_data.find("709_3686")
incMC_709_3686 = DatasetManager.load_inclusive_mc.find("709_3686")

all_data = [data_709_3686]
all_incMC = [incMC_709_3686]

# Decay card: psi(3686) -> Sigma- Sigma+bar
# Sigma- -> n pi-; Sigma+bar -> nbar pi+
# Note: anti-neutron reconstructed via EMC; neutron treated as missing
decay_card = <<~DECAYCARD
    Decay psi(3686)
    1.000  Sigma-  anti-Sigma+                HELAMP 1.0 0.0;
    Enddecay

    Decay Sigma-
    1.000  n0  pi-                             PHSP;
    Enddecay

    Decay anti-Sigma+
    1.000  anti-n0  pi+                        PHSP;
    Enddecay
End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name   = "psi3686_Sigma_Sigmabar_signal"
  config.events         = 500000
  config.decay_card     = decay_card
  config.cross_section  = :default
end

### Event selection ###
alg = Algorithm.new("Psi3686_SigmaSigmabar")
alg.set_header(["Psi3686SigmaAlg/Psi3686SigmaSigmabar.h"])
   .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new

# Charged tracks: exactly 2 (pi+ and pi-), zero net charge
# Anti-neutron reconstructed via EMC (not as a charged track)
event_selection.select_track {
                 cos_theta 0.93
                 Vz   30.0         # |Vz| < 30 cm (looser than standard for Sigma daughters)
                 Vr   10.0         # |Vxy| < 10 cm (looser)
                 nTot "==2"        # 2 charged pions, zero net charge
               }
               # PID: both tracks identified as pions (L(pi) > L(K) and L(pi) > L(p))
               .pid(method: :probability) {
                 prob_cut 0.001
                 identify :pion, against: [:kaon, :proton]
                 npip "==1"
                 npim "==1"
               }
               .assign({:chrgp => :pip, :chrgn => :pim})
               # 4C kinematic fit with neutron as free parameters + Sigma+ mass constraint
               # Anti-neutron angles used (not energy); neutron momentum floated
               .kinematic_fit([:pip, :pim]) {
                 nominal
                 constrain_four_momentum
                 chi2_cut 200
               }

alg.with_decay_card(decay_card).apply(event_selection)

# Anti-neutron reconstruction via EMC showers — inexpressible in DSL
# EMC shower energy > 600 MeV (barrel and endcap)
# Lateral moment (second moment) > 20 to suppress photon background
# Angle to nearest charged track > 10 degrees
# EMC time within [0, 700] ns
# If multiple anti-neutron candidates, select most energetic
alg.note(:antineutron_emc_selection,
  "Anti-neutron reconstructed via EMC showers: E(EMC) > 600 MeV; lateral moment > 20; angle to charged track > 10 deg; EMC time in [0,700] ns; most energetic candidate retained. Applied in ROOT.")

# Kinematic fit with free neutron parameters — inexpressible in DSL
# Neutron three-momentum left as free parameters; anti-neutron polar and azimuthal angles used (energy left free)
# Four-momentum conservation + M(nbar pi+) = M(Sigma+) mass constraint
# chi2 < 50 (optimized by S/sqrt(S+B))
alg.note(:kinematic_fit_free_neutron,
  "Kinematic fit: four-momentum conservation + M(nbar pi+) = M(Sigma+) mass constraint. Neutron three-momentum free parameters; anti-neutron angles used, energy floated. chi2 < 50 (optimized). Applied in ROOT using bes kinematic fit.")

# pi+pi- recoil mass < 2.9 GeV/c^2 to suppress psi(3686) -> pi+pi- J/psi with J/psi -> n nbar
alg.note(:pipi_recoil_mass_veto,
  "Recoil mass of pi+pi- < 2.9 GeV/c^2 to suppress psi(3686) -> pi+pi- J/psi, J/psi -> n nbar. Applied in ROOT.")

# Peaking backgrounds from psi(3686) -> gamma chi_cJ (chi_cJ -> Sigma- Sigma+) and gamma eta_c
alg.note(:peaking_backgrounds,
  "Peaking backgrounds: psi(3686) -> gamma chi_cJ (J=0,1,2), chi_cJ -> Sigma- Sigma+ (562+/-86 events) and psi(3686) -> gamma eta_c, eta_c -> Sigma- Sigma+ (5+/-1 events). Estimated from dedicated MC samples. Applied in ROOT.")

# Off-resonance data at 3.65 GeV to estimate non-psi(3686) background
alg.note(:non_psi3686_background,
  "Non-psi(3686) background estimated from off-resonance data at 3.65 GeV: 92+/-53 events (scaled by luminosity, 1/s^2 and efficiency). Applied in ROOT.")

# Signal extraction via unbinned ML fit to M(n pi-) in [1.15, 1.25] GeV/c^2
# Signal: MC shape convolved with Gaussian
# Peaking BG: fixed from MC exclusive background shapes
# Non-peaking BG: first-order polynomial
alg.note(:signal_extraction,
  "Signal from unbinned ML fit to M(n pi-) in [1.15,1.25] GeV/c^2. Signal: MC shape conv. Gaussian; peaking BG: fixed MC shapes; non-peaking BG: 1st-order polynomial. Applied in ROOT.")

# Angular distribution fit: 1 + alpha_Sigma * cos^2(theta_Sigma) in 10 bins of cos(theta_Sigma)
alg.note(:angular_distribution,
  "Angular parameter alpha_Sigma- measured from efficiency-corrected cos(theta_Sigma-) distribution fit to 1+alpha*cos^2(theta). Applied in ROOT.")

# 448.1M psi(3686) events total
alg.note(:total_psi3686_events,
  "Total of (448.1 +/- 2.9) x 10^6 psi(3686) events used. B(psi(3686) -> Sigma- Sigma+bar) = (2.82 +/- 0.04_stat +/- 0.08_syst) x 10^-4.")

all_datasets = all_data + all_incMC + exMC_signal
root_files = alg.execute_on(all_datasets)