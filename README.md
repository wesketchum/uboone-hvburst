# uboone-hvburst

Setup:
<normal setup of products area>
setup uboonecode v06_23_00 -qe10:prof
setup gallery v1_03_08 -qe10:nu:prof

Make:
make

Run:
./analyzeBurst <list_of_artroot_files.list> <output_file_name.root>

Ntuple produced (as of Feb10):
TNtuple ntp("ntp","Burst Ntuple","run:ev:i_ev:ch:ped:max:max_tick:risetime:falltime:min:min_tick:ref_tick");

Note, i_ev is a per-job event counter, to allow selecting singular events without needing the event number memorized.

See docdb 7071
