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
  TNtuple ntp("ntp","Burst Ntuple","run:ev:i_ev:ch:ped:max:max_tick:risetime:falltime:min:min_tick:ref_tick");
  
  vector<string> filenames;

  //char file_name[256];
  string file_name;
  ifstream input_file(argc[1]);
  while(getline(input_file,file_name))
    filenames.push_back(file_name);
  InputTag rawdigit_tag { "daq" };
  const size_t WINDOW_SIZE = 30;

  size_t ev_counter = 0;
  
  for (gallery::Event ev(filenames) ; !ev.atEnd(); ev.next()) {
    auto t_begin = high_resolution_clock::now();
    
    //to get run and event info, you use this "eventAuxillary()" object.
    cout << "Processing "
	 << "Run " << ev.eventAuxiliary().run() << ", "
	 << "Event " << ev.eventAuxiliary().event() << endl;

    auto const& rawdigit_handle = ev.getValidHandle< vector<raw::RawDigit> >(rawdigit_tag);
    vector<float> pedestal_vec(rawdigit_handle->size());
    vector<short> max_ticks;
    size_t i_ch=0;
    for (auto const& wvfm : *rawdigit_handle){

      /*
      cout << "\tChannel : " << wvfm.Channel();
      auto i_max = std::max_element(wvfm.ADCs().begin(),wvfm.ADCs().end());
      cout << "\t max is " << *i_max << " at " << std::distance(wvfm.ADCs().begin(),i_max) << endl;
      */

      auto i_max = std::max_element(wvfm.ADCs().begin(),wvfm.ADCs().end());
      max_ticks.push_back(std::distance(wvfm.ADCs().begin(),i_max));

      pedestal_vec[i_ch] = std::accumulate(wvfm.ADCs().begin(),wvfm.ADCs().end(),0.0) / (float)(wvfm.ADCs().size());
      //cout << "\tPedestal = " << pedestal_vec[i_ch] << endl;
      ++i_ch;
    }

    std::sort(max_ticks.begin(),max_ticks.end());
    const size_t MAX_TICK = max_ticks.at(max_ticks.size()/2);

    cout << "Using max tick = " << MAX_TICK << endl;
    
    if(MAX_TICK<WINDOW_SIZE || (rawdigit_handle->at(0).ADCs().size()-MAX_TICK)<WINDOW_SIZE){
      cout << "ERROR: MAX_TICK " << MAX_TICK << " too close to edge." << endl;
      return -1;
    }

    i_ch=0;
    for (auto const& wvfm : *rawdigit_handle){

      auto i_max = std::max_element(wvfm.ADCs().begin()+(MAX_TICK-WINDOW_SIZE),
				    wvfm.ADCs().begin()+(MAX_TICK+WINDOW_SIZE));
      auto i_min = std::min_element(wvfm.ADCs().begin()+(MAX_TICK-WINDOW_SIZE),
				    wvfm.ADCs().begin()+(MAX_TICK+WINDOW_SIZE));

      /*
      cout << "\tChannel " << wvfm.Channel() << ":" << endl;
      cout << "\t\t max is " << (*i_max)-pedestal_vec[i_ch] << " at " << std::distance(wvfm.ADCs().begin(),i_max) << endl;
      cout << "\t\t min is " << (*i_min)-pedestal_vec[i_ch] << " at " << std::distance(wvfm.ADCs().begin(),i_min) << endl;
      */

      size_t i_rise = std::distance(wvfm.ADCs().begin(),i_max);
      float thresh = 0.1 * ((float)(*i_max) - pedestal_vec[i_ch]);
      while((wvfm.ADCs()[i_rise]-pedestal_vec[i_ch])>thresh && i_rise>MAX_TICK-WINDOW_SIZE)	--i_rise;

      size_t i_fall = std::distance(wvfm.ADCs().begin(),i_max);
      while((wvfm.ADCs()[i_fall]-pedestal_vec[i_ch])>thresh && i_fall<MAX_TICK+WINDOW_SIZE) ++i_fall;
      
      ntp.Fill(ev.eventAuxiliary().run(),
	       ev.eventAuxiliary().event(),
	       ev_counter,
	       wvfm.Channel(),
	       pedestal_vec[i_ch],
	       *i_max,std::distance(wvfm.ADCs().begin(),i_max),
	       std::distance(wvfm.ADCs().begin(),i_max)-i_rise,i_fall-std::distance(wvfm.ADCs().begin(),i_max),
	       *i_min,std::distance(wvfm.ADCs().begin(),i_min),
	       MAX_TICK);
      ++i_ch;
    }
    
    ++ev_counter;
    
    auto t_end = high_resolution_clock::now();
    duration<double,std::milli> time_total_ms(t_end-t_begin);
    cout << "\tEvent " <<  ev_counter << " took " << time_total_ms.count() << " ms to process." << endl;
  } //end loop over events!


  //and ... write to file!
  f_output.cd();
  f_output.Write();
  f_output.Close();

}