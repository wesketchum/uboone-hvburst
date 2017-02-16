/*************************************************************
 *************************************************************/


//some standard C++ includes
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <algorithm>
#include <chrono>

//some ROOT includes
#include "TInterpreter.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TClonesArray.h"

//"art" includes (canvas, and gallery)
#include "canvas/Utilities/InputTag.h"
#include "gallery/Event.h"
#include "gallery/ValidHandle.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOne.h"

//"larsoft" object includes
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RawData/TriggerData.h"

//our own includes!
//#include "hist_utilities.h"

//convenient for us! let's not bother with art and std namespaces!
using namespace art;
using namespace std;

using namespace std::chrono;

int main(int argv, char** argc) {

  if(argv!=3){
    cout << "ERROR: Usage is analyzeBurst <input_file_list> <output_file>" << endl;
    return -1;
  }
  
  TFile f_output(argc[2],"RECREATE");
  TNtuple ntp("ntp","Burst Ntuple","run:ev:i_ev:ts:ch:ped:max:max_tick:risetime:falltime:min:min_tick:ref_tick");
  TNtuple nt_op("nt_op","OpDet Burst Ntuple","run:ev:i_ev:ts:ch:time_us:ref_tick:rel_time");
  
  vector<string> filenames;

  //char file_name[256];
  string file_name;
  ifstream input_file(argc[1]);
  while(getline(input_file,file_name))
    filenames.push_back(file_name);
  InputTag daq_tag { "daq" };
  InputTag opdet_tag { "pmtreadout","OpdetCosmicHighGain" };
  const size_t WINDOW_SIZE = 30;
  const float  BURST_THRESHOLD = 115.;

  size_t ev_counter = 0;
  
  for (gallery::Event ev(filenames) ; !ev.atEnd(); ev.next()) {
    auto t_begin = high_resolution_clock::now();
    
    //to get run and event info, you use this "eventAuxillary()" object.
    cout << "Processing "
	 << "Run " << ev.eventAuxiliary().run() << ", "
	 << "Event " << ev.eventAuxiliary().event() << ", " 
	 << "Time " << ev.eventAuxiliary().time().timeHigh() << endl;

    auto const& rawdigit_handle = ev.getValidHandle< vector<raw::RawDigit> >(daq_tag);
    vector<float> pedestal_vec(rawdigit_handle->size());
    vector<short> max_ticks;
    size_t i_ch=0;

    const size_t N_ADCs = rawdigit_handle->at(0).ADCs().size();

    
    //fill pedestals
    for (auto const& wvfm : *rawdigit_handle){
      pedestal_vec[i_ch] = std::accumulate(wvfm.ADCs().begin(),wvfm.ADCs().end(),0.0) / (float)(wvfm.ADCs().size());
      //pedestal_vec[i_ch] = wvfm.GetPedestal();
      ++i_ch;
    }
    
    size_t MAX_TICK=0;
    for (size_t iadc=0; iadc<N_ADCs; ++iadc){

      float adc_sum=0; i_ch=0;
      for (auto const& wvfm : *rawdigit_handle){
	adc_sum += ((float)wvfm.ADCs()[iadc] - pedestal_vec[i_ch])/8256.;
	//adc_sum += ((float)wvfm.ADCs()[iadc] - (float)wvfm.GetPedestal())/8256.;
	//std::cout << "\t\tch=" << i_ch << " adcsum is " << adc_sum << std::endl;
	++i_ch;
      }
      if(adc_sum > BURST_THRESHOLD){
	MAX_TICK = iadc;
	break;
      }
    }

    cout << "Using max tick = " << MAX_TICK << endl;
    
    if(MAX_TICK<WINDOW_SIZE || (rawdigit_handle->at(0).ADCs().size()-MAX_TICK)<WINDOW_SIZE){
      cout << "ERROR: MAX_TICK " << MAX_TICK << " too close to edge." << endl;
    }
    else{
      i_ch=0;
      for (auto const& wvfm : *rawdigit_handle){
	
	auto i_max = std::max_element(wvfm.ADCs().begin()+(MAX_TICK-WINDOW_SIZE),
				      wvfm.ADCs().begin()+(MAX_TICK+WINDOW_SIZE));
	auto i_min = std::min_element(wvfm.ADCs().begin()+(MAX_TICK-WINDOW_SIZE),
				    wvfm.ADCs().begin()+(MAX_TICK+WINDOW_SIZE));
	
	size_t i_rise = std::distance(wvfm.ADCs().begin(),i_max);
	float thresh = 0.1 * ((float)(*i_max) - pedestal_vec[i_ch]);
	while((wvfm.ADCs()[i_rise]-pedestal_vec[i_ch])>thresh && i_rise>MAX_TICK-WINDOW_SIZE)	--i_rise;
	
	size_t i_fall = std::distance(wvfm.ADCs().begin(),i_max);
	while((wvfm.ADCs()[i_fall]-pedestal_vec[i_ch])>thresh && i_fall<MAX_TICK+WINDOW_SIZE) ++i_fall;
	
	ntp.Fill(ev.eventAuxiliary().run(),
		 ev.eventAuxiliary().event(),
		 ev_counter,
		 ev.eventAuxiliary().time().timeHigh(),
		 wvfm.Channel(),
		 pedestal_vec[i_ch],
		 *i_max,std::distance(wvfm.ADCs().begin(),i_max),
		 std::distance(wvfm.ADCs().begin(),i_max)-i_rise,i_fall-std::distance(wvfm.ADCs().begin(),i_max),
		 *i_min,std::distance(wvfm.ADCs().begin(),i_min),
		 MAX_TICK);
	++i_ch;
      }
      
      auto const& opdet_handle = ev.getValidHandle<vector<raw::OpDetWaveform>>(opdet_tag);
      auto const& trig_handle = ev.getValidHandle<vector<raw::Trigger>>(daq_tag);
      
      auto const& opdet_vec(*opdet_handle);
      auto trig_time = trig_handle->at(0).TriggerTime();
      
      
      for(auto const& wvfm : opdet_vec){
	nt_op.Fill(ev.eventAuxiliary().run(),
		   ev.eventAuxiliary().event(),
		   ev_counter,
		   ev.eventAuxiliary().time().timeHigh(),
		   wvfm.ChannelNumber(),
		   wvfm.TimeStamp() - trig_time,
		   MAX_TICK,
		   (wvfm.TimeStamp() - trig_time) - ((0.5*(float)(MAX_TICK))-1600.));
      }
    }
    ++ev_counter;
    
    auto t_end = high_resolution_clock::now();
    duration<double,std::milli> time_total_ms(t_end-t_begin);
    cout << "\tEvent " <<  ev_counter << " took " << time_total_ms.count() << " ms to process." << endl;
  } //end loop over events!


  //and ... write to file!
  f_output.Write();
  f_output.Close();

}
