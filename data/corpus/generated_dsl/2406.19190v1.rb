# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset description ###
# Eight energy points of the scan: √s = 4.128, 4.157, 4.178, 4.189, 4.199, 4.209, 4.219, 4.226 GeV
scan_sample_names = %w[705_4130 705_4160 703_4180 703_4190 703_4200 703_4210 703_4220 703_4230]

data_points = scan_sample_names.map { |name| DatasetManager.real_data.find(name) }    # real data at each scan point
incMCs      = scan_sample_names.map { |name| DatasetManager.inclusive_mc.find(name) } # matching inclusive MC

# Decay card for the signal process
# ψ(4260) → Ds*− Ds+, Ds*− → γ Ds−, Ds+ → K0S e+ νe (ISGW2), K0S → π+π−, tag Ds− → K+K−π−
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.000 D_s*- D_s+ PHSP;
    Enddecay

    Decay D_s*-
    1.000 gamma D_s- PHSP;
    Enddecay

    Decay D_s+
    1.000 K_S0 e+ nu_e ISGW2;
    Enddecay

    Decay D_s-
    1.000 K+ K- pi- PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# 500k-event exclusive signal MC generated for every scan energy point
exMCs_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_ds_k0s_e_nu"
  config.events        = 500_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "DsToK0SeNu"
ds_k0s_e_nu = TagAnalysis.new(alg_name)   # tag-based analysis (uses pre-stored DTag candidates)
ds_k0s_e_nu.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 4.178]})  # central value of the per-run measured E_cms
            .with_decay_card(decay_card_signal)

# Tag side: single tag Ds− (charm = −1) reconstructed in 14 hadronic modes
ds_k0s_e_nu.tag_side(:Ds) do |t|
  t.modes :DstoKKPi,                          # K+K−π−
          :DstoKKPiPi0,                       # K+K−π−π0
          :DstoPiPiPi,                        # π+π−π−
          :DstoKsK,                           # K0S K−
          :DstoKsKPi0,                        # K0S K−π0
          :DstoKPiPi,                         # K−π+π−
          :DstoKsKsPi,                        # K0S K0S π−
          :DstoKsKPiPi,                       # K0S K−π+π−
          :DstoKsKPiPiPi0,                    # K0S K−π+π−π0
          :DstoEtaGamGamPi,                   # η(γγ) π−
          :DstoEtaPrimePiPiEtaRho,            # η′(π+π−η) ρ−
          :DstoEtaPrimeEtaGamGamPiPiPi,       # η′(ηγγπ+π−) π−
          :DstoEtaPrimeGamRhoPi,              # η′(γρ) π−
          :DstoEtaGamGamRho                   # η(γγ) ρ−
  t.charm -1
end

# Signal side: Ds+ → K0S e+ νe, K0S → π+π−; exactly one π+, one π−, one e+, net charge +1, zero photons
ds_k0s_e_nu.signal_side do |s|
  s.charged(ep: 1, pip: 1, pim: 1)  # one positron + the two pions of the K0S
  s.photons 0                       # no photon participates on the signal side
  s.require_charge 1                # net charge of the signal side
  s.missing :nu_e                   # the undetected neutrino (massless)
end

# 4C kinematic fit with an additional K0S mass constraint
ds_k0s_e_nu.fit do |f|
  f.constrain_four_momentum                                      # 4C: constraint to the measured CMS 4-vector
  f.invariant_mass_of(:pip, :pim).constrain_to_nominal_mass_of(:K_S0)  # K0S mass constraint
  f.chi2_cut 200                                                 # χ² < 200
end

# BOSS-side procedures without a dedicated DSL method
ds_k0s_e_nu
  .note(:transition_gamma_hypothesis, "Two transition-γ hypotheses are tested and the smaller-χ² one is retained:
    (a) γ + ST Ds− → Ds*− and (b) γ + semileptonic Ds+ → Ds*+; only the combination with the smaller χ² (< 100) is kept.")
  .note(:positron_id, "Positron identification combines dE/dx, TOF and EMC likelihoods:
    CLe' > 0.001 and CLe'/(CLe'+CLπ'+CLK') > 0.8.")
  .note(:k0s_selection, "K0S candidates are required to satisfy |M(π+π−) − m_K0S| < 12 MeV/c² and a decay length
    larger than twice the vertex resolution.")
  .note(:unused_shower_energy, "The largest energy among showers not used in the reconstruction must be < 0.2 GeV.")
  .note(:background_veto, "Events with M(K0S e+) < 1.78 GeV/c² are vetoed to suppress the Ds+ → K0S K+ background.")
  .note(:pion_momentum, "Pions not originating from K0S, η or η′ are required to have momentum > 0.1 GeV/c.")

ds_k0s_e_nu.apply   # render the tag specification (takes no Selection argument)
root_files = ds_k0s_e_nu.execute_on(data_points + incMCs + exMCs_signal)