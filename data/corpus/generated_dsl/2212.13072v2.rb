# frozen_string_literal: true
### Dataset description ###
data_4180  = DatasetManager.real_data.find("703_4180")       # 4.178 GeV real data (BOSS 703, 3.19 fb^-1)
incMC_4180 = DatasetManager.inclusive_mc.find("703_4180")    # matching inclusive MC sample

# Decay card for the K_S^0 K^- tag mode (K_S^0 -> pi+ pi-) together with the signal D_s^+ -> pi+ pi+ pi- X
decay_card_tag_ksk = <<~DECAYCARD
    Decay psi(4260)
    1.0000  D_s+  D_s-   PHSP;
    Enddecay

    Decay D_s-
    1.0000  K_S0  K-     PHSP;
    Enddecay

    Decay K_S0
    1.0000  pi+  pi-     PHSP;
    Enddecay

    Decay D_s+
    1.0000  pi+  pi+  pi-   PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for the K^- K^+ pi^- tag mode together with the signal D_s^+ -> pi+ pi+ pi- X
decay_card_tag_kkpi = <<~DECAYCARD
    Decay psi(4260)
    1.0000  D_s+  D_s-   PHSP;
    Enddecay

    Decay D_s-
    1.0000  K-  K+  pi-   PHSP;
    Enddecay

    Decay D_s+
    1.0000  pi+  pi+  pi-   PHSP;
    Enddecay

    End
DECAYCARD

# 500k exclusive MC events for each of the two tag modes
exMC_tag_ksk = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_4180_dstag_ksk"
  config.related_dataset = data_4180
  config.events          = 500_000
  config.decay_card      = decay_card_tag_ksk
  config.cross_section   = :default
end

exMC_tag_kkpi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_4180_dstag_kkpi"
  config.related_dataset = data_4180
  config.events          = 500_000
  config.decay_card      = decay_card_tag_kkpi
  config.cross_section   = :default
end

### Event selection (BOSS) — tag-based double-tag measurement ###
alg_name   = "DsTag3PiX"
ds_tag_alg = TagAnalysis.new(alg_name)
ds_tag_alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
          .set_constant({ "ECMS" => [:double, 4.178] })
          .with_decay_card(decay_card_tag_ksk)

# Tag side: D_s^- (charm -1) reconstructed in either K_S^0 K^- or K^- K^+ pi^-
ds_tag_alg.tag_side(:Ds) do |t|
  t.modes :DstoKsK, :DstoKKPi     # the two tag modes
  t.charm -1                      # pin the tagged side to D_s^-
  t.window :mBC, min: 2.05        # user-requested M_BC > 2.05 GeV/c^2
end

# Signal side: two pi+ and one pi- (net charge +1) plus one undetected inclusive X
ds_tag_alg.signal_side do |s|
  s.charged(pip: 2, pim: 1)
  s.require_charge 1
  s.missing :X0, mass: nil        # inclusive undetected state X (massless form)
end

# 4C kinematic fit (four-momentum constraint) with chi^2 < 200
ds_tag_alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

# BOSS-side procedures that cannot be expressed in the tag DSL
ds_tag_alg
  .note(:background_veto, "signal-side K_S0 veto: reject any pi+ pi- pair with invariant mass in [0.485, 0.510] GeV/c^2")
  .note(:electron_veto, "signal-side electron veto: require E/p < 0.8 for every signal charged track")
  .note(:best_tag_selection, "choose the best D_s^- tag candidate by minimum |deltaE|")
  .note(:mbc_window, "the M_BC > 2.05 GeV/c^2 requirement nominally applies only to the K_S^0 K^- tag mode; the single tag-side window declared here is applied to both tag modes")

ds_tag_alg.apply    # no Selection argument

# Execute on real data, inclusive MC and the two exclusive tag-mode MC samples
root_files = ds_tag_alg.execute_on([data_4180, incMC_4180, exMC_tag_ksk, exMC_tag_kkpi])