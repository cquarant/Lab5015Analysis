#!/usr/bin/env python3
"""One bounded raw-data pass: inspect availability of every REF channel pair."""
import argparse
import csv
import json
from pathlib import Path
import ROOT

CODE = r'''
#include <TChain.h>
#include <TH1D.h>
#include <TFile.h>
#include <fstream>
#include <vector>
#include <map>
void refOccupancy(const std::vector<std::string>& files, const std::vector<double>& factors,
                  const std::string& output, long long expected) {
 TChain t("data"); for(auto& f:files) t.Add(f.c_str());
 float ov, step; int idx[256]; std::vector<float>* en=nullptr; std::vector<float>* tot=nullptr;
 t.SetBranchStatus("*",0);
 for(auto b:{"step1","step2","channelIdx","energy","tot"}) t.SetBranchStatus(b,1);
 t.SetBranchAddress("step1",&ov); t.SetBranchAddress("step2",&step);
 t.SetBranchAddress("channelIdx",idx); t.SetBranchAddress("energy",&en); t.SetBranchAddress("tot",&tot);
 const int m[]={4,29,1,25,0,24,3,30,2,31,6,28,7,27,9,22,5,26,11,21,8,23,12,20,10,19,14,18,15,16,13,17};
 TFile f((output+".root").c_str(),"RECREATE");
 std::map<int,TH1D*> h; std::map<int,long long> pairs, selected, dut7;
 long long total=t.GetEntries();
 if(total!=expected) throw std::runtime_error("Input count differs from spill manifest");
 for(long long i=0;i<total;++i) {
  if(t.GetEntry(i)<=0) throw std::runtime_error("Cannot read raw entry");
  if(ov!=3.f) continue; int th=int(step/10000)-1; int active=0;
  for(int b=0;b<16;++b) {
   int l=idx[128+m[2*b]],r=idx[128+m[2*b+1]]; if(l<0||r<0) continue;
   if((*tot)[l] > -10000 && (*tot)[l] < 100000 && (*tot)[r] > -10000 && (*tot)[r] < 100000 &&
      ((*en)[l]+(*en)[r])>0) ++active;
  }
  for(int b=0;b<16;++b) {
   int k=100*th+b,l=idx[128+m[2*b]],r=idx[128+m[2*b+1]]; if(l<0||r<0) continue;
   ++pairs[k]; if(active>3) continue; ++selected[k];
   if(!h.count(k)) h[k]=new TH1D(Form("energy_ref%02d_th%02d",b,th),"Calibrated REF mean energy",512,0,1024);
   h[k]->Fill(0.5*((*en)[l]*factors[2*b]+(*en)[r]*factors[2*b+1]));
   if(idx[224+m[14]]>=0 && idx[224+m[15]]>=0) ++dut7[k];
  }
 }
 std::ofstream out(output+".csv"); out<<"ref,threshold,pair_present,after_multiplicity,dut7_present,raw_events\n";
 for(int th:{3,5,7,10,15,20}) for(int b=0;b<16;++b) { int k=100*th+b;
  out<<b<<","<<th<<","<<pairs[k]<<","<<selected[k]<<","<<dut7[k]<<","<<total<<"\n";
 }
 f.Write(); f.Close();
}
'''

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--baseline',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(); a.output.mkdir(parents=True,exist_ok=True)
    ROOT.gROOT.SetBatch(True)
    if not ROOT.gInterpreter.Declare(CODE): raise RuntimeError('Cannot compile raw audit')
    spills=list(csv.DictReader((a.baseline/'input_spills.csv').open()))
    for run in (4194,4257,4281,4291):
        dst=a.output/str(run)
        expected=sum(int(r['entries']) for r in spills if int(r['run'])==run)
        if dst.with_suffix('.csv').exists() and dst.with_suffix('.root').exists():
            old=list(csv.DictReader(dst.with_suffix('.csv').open()))
            if len(old)==96 and all(int(r['raw_events'])==expected for r in old): continue
            raise RuntimeError('Incomplete existing audit; use a new output directory')
        calib=list(csv.DictReader((a.baseline/str(run)/'TOFHIR_REF_calibration_factors.csv').open()))
        factors={(int(r['bar']),r['side']):float(r['calib']) for r in calib if float(r['vov'])==3 and int(r['th'])==10}
        ROOT.refOccupancy([r['reco_path'] for r in spills if int(r['run'])==run],
                          [factors[b,s] for b in range(16) for s in ('L','R')], str(dst), expected)
        print(run,'complete',flush=True)

if __name__=='__main__': main()
