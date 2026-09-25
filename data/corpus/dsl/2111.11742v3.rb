DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

psi3770_data  = DatasetManager.real_data.find("712_3773")
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

# Decay card: e+ e- -> Lambda anti-Lambda at psi(3770)
# One-photon + psi(3770) exchange; use KKMC with psi(4260) top mother
decay_card_lalabar = <<~DECAYCARD
    Decay psi(4260)
    1.000 Lambda0 anti-Lambda0 PHSP;
    Enddecay
    Decay Lambda0
    1.000 p+ pi- HypWK;
    Enddecay
    Decay anti-Lambda0
    1.000 anti-p- pi+ HypWK;
    Enddecay
    End
DECAYCARD

exMC_lalabar = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "signal_lalabar_3773"
  config.related_dataset = psi3770_data
  config.events          = 200_000_000
  config.decay_card      = decay_card_lalabar
  config.cross_section   = :default
end

alg_lalabar = Algorithm.new("LambdaLambdaBar")
alg_lalabar.set_header(["LambdaLambdaBarAlg/LambdaLambdaBar.h"])
            .set_constant({"ECMS" => [:double, 3.773]})
            .note(:decay_length_cut,
              "Lambda decay length > 0 required to suppress non-Lambda background " \
              "(negative decay lengths are from detector resolution).")
            .note(:mass_window_cut,
              "M(p pi-) and M(p-bar pi+) required within 5 MeV/c^2 of known Lambda mass. " \
              "Signal region determined by figure-of-merit S/sqrt(S+B) from MC. " \
              "Background estimated by corner method from 4 sideband regions.")

event_sel = Selection.new
  .select_track do
    cos_theta   0.93
    Vr          1.0
    Vz          10.0
    nChrp       "==2"
    nChrn       "==2"
    nNet        "==0"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprp ">=1"
    nprm ">=1"
  end
  .remove([:prp <= :chrgp, :prm <= :chrgn])
  .assign({chrgp: :pip, chrgn: :pim})
  .secondary_vertex_fit([:prp, :pim]) do
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  .secondary_vertex_fit([:prm, :pip]) do
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  .kinematic_fit([:Lambda, :Lambda_bar]) do
    nominal
    constrain_four_momentum
    chi2_cut 200
  end

alg_lalabar.with_decay_card(decay_card_lalabar).apply(event_sel)
alg_lalabar.execute_on([psi3770_data, psi3770_incMC, exMC_lalabar])