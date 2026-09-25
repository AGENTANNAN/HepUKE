### Dataset description ###
# Cross-section scan: 48 CM energies between 3.51 and 4.95 GeV (44.2 fb^-1 total).
datasets_scan = [
  DatasetManager.real_data.find("712_3773"),
  DatasetManager.real_data.find("703_4180"),
  DatasetManager.real_data.find("703_4190"),
  DatasetManager.real_data.find("703_4200"),
  DatasetManager.real_data.find("703_4210"),
  DatasetManager.real_data.find("703_4220"),
  DatasetManager.real_data.find("703_4230"),
  DatasetManager.real_data.find("703_4237"),
  DatasetManager.real_data.find("703_4246"),
  DatasetManager.real_data.find("703_4260"),
  DatasetManager.real_data.find("703_4270"),
  DatasetManager.real_data.find("703_4280"),
  DatasetManager.real_data.find("703_4360"),
  DatasetManager.real_data.find("703_4420"),
  DatasetManager.real_data.find("703_4600"),
  DatasetManager.real_data.find("706_4610"),
  DatasetManager.real_data.find("706_4620"),
  DatasetManager.real_data.find("706_4640"),
  DatasetManager.real_data.find("706_4660"),
  DatasetManager.real_data.find("706_4680"),
  DatasetManager.real_data.find("706_4700"),
  DatasetManager.real_data.find("707_4740"),
  DatasetManager.real_data.find("707_4750"),
  DatasetManager.real_data.find("707_4780"),
  DatasetManager.real_data.find("707_4840"),
  DatasetManager.real_data.find("707_4914"),
  DatasetManager.real_data.find("707_4946"),
]

# ConExc-style continuum decay card for Xi(1530)0 anti-Xi0 production
signal_decay_card = <<~DECAYCARD
  Decay vpho
  1.000 Xi*0 anti-Xi0                      PHSP;
  Enddecay

  Decay Xi*0
  1.000 Xi- pi+                            PHSP;
  Enddecay

  Decay Xi-
  1.000 Lambda0 pi-                        PHSP;
  Enddecay

  Decay Lambda0
  1.000 p+ pi-                             PHSP;
  Enddecay

  Decay anti-Xi0
  1.000 anti-Lambda0 pi0                   PHSP;
  Enddecay

  Decay anti-Lambda0
  1.000 anti-p- pi+                        PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma                        PHSP;
  Enddecay

  End
DECAYCARD

# Exclusive MC per energy point
exMC_signal_list = DatasetManager.create_exclusive_mc_for(datasets_scan) do |c|
  c.sample_name    = "Xi1530_Xi0_signal"
  c.events         = 100_000
  c.decay_card     = signal_decay_card
  c.cross_section  = :default
end
exMC_signal_list.each { |m| m.save_to_config(format: :yaml, file_path: 'temp_for_test') }

### Event selection ###
alg_name = "EEToXi1530Xi0"
alg = Algorithm.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({ "ECMS" => [:double, 3.773] })
   .set_alias({ "std::vector<double>" => "Vdouble" })

sel = Selection.new
sel.select_track {
      cos_theta 0.93
      nChrp ">=2"
      nChrn ">=2"
    }
   .pid(method: :probability) {
      prob_cut 0.001
      identify :proton, against: [:kaon, :pion]
      identify :pion,   against: [:kaon, :proton]
      nprp ">=1"
      nprm ">=1"
      npip ">=1"
      npim ">=2"
    }
   # Lambda reconstruction: p pi-
   .secondary_vertex_fit([:prp, :pim]) {
      build_virtual_particle(:Lambda).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
   }
   # Xi- reconstruction: Lambda pi-
   .secondary_vertex_fit([:Lambda, :pim]) {
      build_virtual_particle(:"Xi-").by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
   }
   # Xi(1530)^0 reconstruction: Xi- pi+  (kinematic fit as bookkeeping; 4C on visible + missing recoil)
   .kinematic_fit([:"Xi-", :pip]) {
      nominal
      constrain_four_momentum
      chi2_cut 200
   }

alg.note(:lambda_mass_window,
         "|M(p pi-) - m_Lambda| < 7 MeV/c^2 and secondary vertex fit chi2 < 500, decay length > 0.")
   .note(:xi_mass_window,
         "|M(Lambda pi-) - m_Xi-| < 6.5 MeV/c^2, secondary vertex fit, decay length > 0, " \
         "best candidate by closest mass to nominal Xi-.")
   .note(:xi1530_selection,
         "Xi(1530)^0 reconstructed from Xi- pi+ by minimising |M(Xi- pi+) - m_Xi(1530)|; " \
         "signal window |M(Xi- pi+) - m_Xi(1530)^0| < 13 MeV/c^2 applied at ROOT stage.")
   .note(:anti_xi0_recoil_window,
         "anti-Xi^0 tagged via recoil against Xi(1530)^0: |M_recoil(Xi- pi+) - m_anti-Xi0| < 25 MeV/c^2 " \
         "at ROOT stage. Signal region [1.277, 1.352] GeV/c^2, sideband [1.440, 1.490] and " \
         "[1.165, 1.190] GeV/c^2.")
   .note(:xi1530_sideband,
         "Xi(1530)^0 sideband [1.602, 1.658] GeV/c^2 used for counting-method background subtraction " \
         "with scale factor R2 = 0.52.")
   .note(:single_baryon_tag,
         "Single-baryon tagging: reconstruct Xi(1530)^0 side only; anti-Xi^0 partner extracted from " \
         "recoil mass spectrum.")

alg.with_decay_card(signal_decay_card).apply(sel)
alg.execute_on(datasets_scan + exMC_signal_list)
