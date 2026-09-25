### Dataset description ###
# Seven energy points from 4.600 to 4.699 GeV (4.5 fb^-1)
# 4600 lives under BOSS 703; 4610-4700 under BOSS 706
sample_names = ["703_4600", "706_4610", "706_4620", "706_4640",
                "706_4660", "706_4680", "706_4700"]

# Real data and inclusive MC at each of the seven points
data_points  = sample_names.map { |n| DatasetManager.real_data.find(n) }
incmc_points = sample_names.map { |n| DatasetManager.inclusive_mc.find(n) }

# Decay card for the signal e+e- -> Lambda_c+ Lambda_c-,
# Lambda_c+ -> Sigma- K+ pi+, Sigma- -> n pi- (neutron undetected)
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0000  Lambda_c+   anti-Lambda_c-   PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000  Sigma-  K+  pi+   PHSP;
    Enddecay

    Decay Sigma-
    1.0000  n0  pi-   PHSP;
    Enddecay

    End
DECAYCARD

# 100k-event exclusive MC for the signal, generated at every energy point
exMCs = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "Lc_SigmaKPi_signal"
  config.events        = 100_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) — tag-based Lambda_c+ Lambda_c- analysis ###
alg_name = "LcSigmaKPi"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 4.600]})   # beam energy is read per-run from the DB
   .with_decay_card(decay_card_signal)

# Tag side: Lambda_c- reconstructed from the pre-stored DTag candidates.
# Mode names use the full DTagAlg channel-name symbols (Lambda_c+ convention);
# charm -1 pins the tagged side to Lambda_c-.
alg.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoKPiP,        # p K- pi+      (desc: p K+ pi-)
          :LambdacPtoKsP,         # p K_S0        (desc: p K_S0)
          :LambdacPtoKPiPi0P,     # p K- pi+ pi0  (desc: p K+ pi- pi0)
          :LambdacPtoKsPi0P,      # p K_S0 pi0    (desc: p K_S0 pi0)
          :LambdacPtoKsPiPiP,     # p K_S0 pi+ pi-(desc: p K_S0 pi+ pi-)
          :LambdacPtoLambdaPi,    # Lambda pi+    (desc: Lambda pi-)
          :LambdacPtoLambdaPiPi0, # Lambda pi+ pi0
          :LambdacPtoLambdaPiPiPi,# Lambda pi+ pi- pi+
          :LambdacPtoSigma0Pi     # Sigma0 pi+    (desc: Sigma0 pi-)
  t.charm -1
end

# Signal side: the tracks/showers the tag did not use.
# Lambda_c+ -> Sigma- K+ pi+, Sigma- -> n pi- (neutron missing):
#   charged multiset  K+ , pi+ , pi-  (total charge +1)
#   neutron declared as a massive missing particle
alg.signal_side do |s|
  s.charged(kp: 1, pip: 1, pim: 1)   # exact multiplicity -> no extra charged tracks allowed
  s.require_charge 1                 # K+ + pi+ + pi- = +1
  s.missing :n                       # neutron treated as missing (massive)
end

# Kinematic fit: 4-momentum conservation + Lambda_c mass constraint (chi2 < 200).
# The neutron participates in the mass constraint (all Lambda_c+ daughters).
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:kp, :pip, :pim, :n).constrain_to_nominal_mass_of(:"Lambda_c+")
  f.chi2_cut 200
end

# Procedures that have no DSL expression on the tag-based surface
alg.note(:signal_track_quality,
         'signal-side K+ and pi+ required |cos theta|<0.93 with |Vz|<10 cm, |Vr|<1 cm; the pi- from Sigma- is allowed |Vz|<20 cm; no additional charged track may pass the track quality')
   .note(:pid_correction_method,
         'signal-side K+/pi+ identified with combined dE/dx and TOF PID probability')
   .note(:signal_vertex_fit,
         'a common vertex fit on the signal K+ and pi+ is applied before the 4C kinematic fit')
   .note(:background_veto,
         'recoil mass against Lambda_c- within (2.275,2.310) GeV/c^2; recoil mass against H+ (Lambda_c+ -> Sigma+ K+ pi-) outside (1.15,1.24) GeV/c^2; M(pi+pi-) outside (0.48,0.52) GeV/c^2 (K_S0 veto); M_rec(H-) > 1.15 GeV/c^2')

alg.apply
alg.execute_on(data_points + incmc_points + exMCs)