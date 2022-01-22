#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <map>

#include <stdio.h>
#include <stdlib.h>

#include <TCanvas.h>
#include <TGraph.h>
#include <TF1.h>
#include <TH1D.h>

#define LABLECALISTART "<<CALIBRATION>>"
#define LABLECALIEND "<<DATA>>"
#define LABELDATASTART "<<DATA>>"
#define LABELDATAEND "<<END>>"

#define EFFENERGY1 1.85
#define EFFENERGY2 15.

#define DEBUG

using namespace::std;

void ReadXRayData(map<string, map<double, double>> &m_p);
void PreAnaXRayData_1(map<string, map<double, double>> &m_p);
void PreAnaXRayData_2(map<string, map<double, double>> &m_p);
void ReadEffData(double par1[10], double par2[10], double par3[10]);
void ReadMcaData(char filename[1024], double par[2], vector<double> &v_x, vector<double> &v_y);
void GetPeakInfo(double par[2], vector<double> &v_x, vector<double> &v_y, map<double, double> &m_p);
void CaliDetectorEff(double par1[10], double par2[10], double par3[10], map<double, double> &m_p, map<double, double> &m_pp);
void CaliAttenEff(map<double, double> &m_p, map<double, double> &m_pp);
void GetElementShow(vector<string> &v_e, map<string, map<double, double>> &m_p, map<double, double> &m_pp, map<string, map<double, double>> &m_ppp, map<string, map<double, double>> &m_q);
void CheckElementShow(map<string, map<double, double>> &m_p, map<string, map<double, double>> &m_pp, map<string, map<double, double>> &m_ppp, map<string, map<double, double>> &m_q);
void GetElementPercent(map<string, map<double, double>> &m_p, map<string, map<double, double>> &m_pp, map<string, double> &m_q);

//
void analysis()
{
  gROOT->SetBatch(0);

  map<string, map<double, double>> map_xray_data;
  ReadXRayData(map_xray_data);
  PreAnaXRayData_1(map_xray_data);
  map<string, map<double, double>>::iterator it1;
  map<double, double>::iterator it2;

  for(it1=map_xray_data.begin();it1!=map_xray_data.end();++it1){
    cout << it1->first << " => \n";
    for(it2=it1->second.begin();it2!=it1->second.end();++it2){
      cout << it2->first << "  " << it2->second << "  "; 
    }
    cout << "\n";
  }

  //need to choose here
  //1
  //vector<string> v_element = {"20Ca", "26Fe", "82Pb", "25Mn", "24Cr", "29Cu", "30Zn", "22Ti"};
  //2
  /*
  v_element.clear();
  for(it1=map_xray_data.begin();it1!=map_xray_data.end();++it1){
    v_element.push_back(it1->first);
  }
  */

  /*
  double par1_eff[10];
  double par2_eff[10];
  double par3_eff[10];
  ReadEffData(par1_eff, par2_eff, par3_eff);
  double par_cali[2];
  vector<double> v_x;
  vector<double> v_y;
  char file_name[1024];
  sprintf(file_name, "../spectrum/CaTi-10KV-0.02mA.mca");
  ReadMcaData(file_name, par_cali, v_x, v_y);
  map<double, double> m_peak;
  map<double, double> m_peak_dector_eff;
  map<double, double> m_peak_dector_eff_atten_eff;
  GetPeakInfo(par_cali, v_x, v_y, m_peak);
  CaliDetectorEff(par1_eff, par2_eff, par3_eff, m_peak, m_peak_dector_eff);
  CaliAttenEff(m_peak_dector_eff, m_peak_dector_eff_atten_eff);

  for(it2=m_peak.begin();it2!=m_peak.end();++it2){
    cout << ",,, raw " << it2->first << " ==>  " << it2->second << endl; 
  }
  for(it2=m_peak_dector_eff.begin();it2!=m_peak_dector_eff.end();++it2){
    cout << ",,, after detector eff " << it2->first << " ==>  " << it2->second << endl; 
  }
  for(it2=m_peak_dector_eff_atten_eff.begin();it2!=m_peak_dector_eff_atten_eff.end();++it2){
    cout << ",,, after detector&atten eff " << it2->first << " ==>  " << it2->second << endl; 
  }

  map<string, map<double, double>> m_element_show;
  map<string, map<double, double>> m_element_show_mca_toi_associate;
  GetElementShow(v_element, map_xray_data, m_peak_dector_eff_atten_eff, m_element_show, m_element_show_mca_toi_associate);

  map<string, map<double, double>> m_element_show_check;
  CheckElementShow(map_xray_data, m_element_show, m_element_show_mca_toi_associate, m_element_show_check);
  map<string, double> m_percent_result;
  GetElementPercent(map_xray_data, m_element_show_check, m_percent_result);
  */
}

//////
void ReadXRayData(map<string, map<double, double>> &m_p)
{
  //m_p: toi x-ray info
  ifstream file_in("../X-Ray.dat");
  if(!file_in){
    std::cout << "can not open X-Ray.dat file" << std::endl;
    return;
  }

  std::string line;
  std::string element;
  double energy, intensity;
  while(1){
    getline(file_in, line);
    std::istringstream iss(line);
    if(!file_in.good()) break;
    iss >> element;
    m_p.insert(pair<string, map<double, double>>(element, map<double,double>()));
    while(iss >> energy >> intensity){
      //cout << energy << " " << intensity << endl;
      m_p[element].insert(pair<double, double>(energy, intensity));
    }
  }

  file_in.close();

  return;
}

//
//only keep the strongest peak in the x-ray data
//m_p: toi x-ray info
void PreAnaXRayData_1(map<string, map<double, double>> &m_p)
{
  map<string, map<double, double>>::iterator it;

  map<double, double> m_pp;//x-ray intensity & x-ray energy
  for(it=m_p.begin();it!=m_p.end();++it){
    m_pp.clear();
    int size = it->second.size();
    if(size > 1){
      for(map<double, double>::iterator itt=it->second.begin();itt!=it->second.end();++itt){
        m_pp.insert(pair<double, double>(itt->second, itt->first));
      }
#ifdef DEBUG
      for(map<double, double>::iterator itt=m_pp.begin();itt!=m_pp.end();++itt){
        cout << it->first << " ==> " << itt->first << " ==> " << itt->second << endl;
      }
#endif
      (it->second).clear();
      map<double, double>::reverse_iterator rit;
      rit=m_pp.rbegin();
      (it->second).insert(pair<double, double>(rit->second, rit->first));

    }else{
      continue;
    }
  }
}

//
//merge peaks energy difference smaller than 0.3keV in the x-ray data
void PreAnaXRayData_2(map<string, map<double, double>> &m_p)
{
  //m_p: toi x-ray info
  map<string, map<double, double>>::iterator it;

  map<double, int> m_diff;//need to be sorted
  for(it=m_p.begin();it!=m_p.end();++it){
    while(1){
      int size = it->second.size();
      //cout << "size  " << size << endl;
      if(size > 1){
        //cout << it->first << endl;

        vector<double> v_energy;
        //
        for(map<double, double>::iterator itt=it->second.begin();itt!=it->second.end();++itt){
          v_energy.push_back(itt->first);
        }
        sort(v_energy.begin(), v_energy.end());
        /*
        cout << "energy ";
        for(double d:v_energy){
          cout << d << " ";
        }
        cout << endl;
        */

        map<double, int> m_diff;
        for(int i=1;i<v_energy.size();i++){
          m_diff.insert(pair<double, int>(v_energy[i]-v_energy[i-1], i));
        }

        map<double, int>::iterator itu = m_diff.begin();
        //cout << " ,, " << itu->first << " " << itu->second << endl;
        if(itu->first<0.3){
          //cout << itu->first << " ==> " << itu->second << " !! " << v_energy[itu->second] << "  " << v_energy[itu->second-1] << endl; 
          //cout << " !! " << v_energy[itu->second] << " " << it->second[v_energy[itu->second]] << "  " << v_energy[itu->second-1] << " " << it->second[v_energy[itu->second-1]] << endl;
          double ea = v_energy[itu->second];
          double ia = it->second[v_energy[itu->second]];
          double eb = v_energy[itu->second-1];
          double ib = it->second[v_energy[itu->second-1]];
          double ee = ea*(ia/(ia+ib))+eb*(ib/(ia+ib));
          double ii = ia+ib;
          it->second.erase(ea);
          it->second.erase(eb);
          it->second.insert(pair<double, double>(ee, ii));
        }
        else  break;
        /*
        for(map<double, double>::iterator itt=it->second.begin();itt!=it->second.end();++itt){
          cout << " ?? " << itt->first << "  " << itt->second << endl;
        }
        */
      }
      else  break;
    }//while
  }
}

//
void ReadEffData(double par1[10], double par2[10], double par3[10])
{
  //par1, par2, par3: curve of 3 cubric function for detector efficiency
  //read effdata
  double par_temp1, par_temp2, par_temp3;
  ifstream ifs_eff;
  ifs_eff.open("eff.par");
  if(!ifs_eff){
    cout << "cannot open eff.par file !" << endl;
    return;
  }
  int kk = 0;
  while(1){
    ifs_eff >> par_temp1 >> par_temp2 >> par_temp3;
    if(!ifs_eff.good()) break;
    par1[kk] = par_temp1;
    par2[kk] = par_temp2;
    par3[kk] = par_temp3;
    kk++;
  }
  ifs_eff.close();
  /*
  for(int i=0;i<10;i++){
    cout << par1[i] << " " << par2[i] << " " << par3[i] << endl;
  }
  */
}

//
void ReadMcaData(char filename[1024], double par[2], vector<double> &x, vector<double> &y)
{
  //filename: mca file used to analysis
  //par: calibration par
  //x:channel info
  //y:counts info
  ifstream fi;
  fi.open(filename);
  if(!fi){
    cout << "cannot open mca file ..." << endl;
    return -1;
  }

  string temp;

  vector<double> cali_ch;
  vector<double> cali_e;
  double cali_ch_ = 0.;
  double cali_e_ = 0.;

  while(1){
    temp.clear();
    getline(fi, temp);
    temp.erase(temp.find_last_not_of("\n\r") + 1);
    if (!fi.good())   break;
    //cout << temp << endl;

    //read cali points
    if(temp.compare(LABLECALISTART) == 0){
      getline(fi, temp);
      temp.erase(temp.find_last_not_of("\n\r") + 1);
      while(1){
        getline(fi, temp);
        temp.erase(temp.find_last_not_of("\n\r") + 1);

        if(temp.compare(LABLECALIEND) == 0){
          break;
        }

        stringstream ss;
        ss.str("");
        ss << temp;
        ss >> cali_ch_ >> cali_e_;
        cali_ch.push_back(cali_ch_);
        cali_e.push_back(cali_e_);
      }
    }

    //read data
    if(temp.compare(LABELDATASTART) == 0){
      cout << temp << endl;
      double xx = 1;
      while(1){
        getline(fi, temp);
        temp.erase(temp.find_last_not_of("\n\r") + 1);
        if(temp.compare(LABELDATAEND) == 0){
          //cout << temp << endl;
          break;
        }
        x.push_back(xx);
        y.push_back((double)atoi(temp.c_str()));
        xx += 1.;
      }
    }
  }

  //
  //cout << cali_ch.size() << endl;
  //cout << x.size() << endl;

  TGraph* gcali = new TGraph();
  for(int i=0;i<cali_ch.size();i++) {
    gcali->SetPoint(i, cali_ch[i], cali_e[i]);
  }

  TCanvas *cc = new TCanvas("cc", "cc", 0, 0, 400, 360);
  cc->cd();
  gcali->GetXaxis()->SetTitle("Channel");
  gcali->GetYaxis()->SetTitle("Energy[keV]");
  gcali->Draw("AP*");

  TF1 *tf1 = new TF1("tf1", "[0]+[1]*x");
  gcali->Fit("tf1", "Q");
  //cout << tf1->GetParameter(0) << endl;
  //cout << tf1->GetParameter(1) << endl;
  par[0] = tf1->GetParameter(0);
  par[1] = tf1->GetParameter(1);

  fi.close();
  delete gcali;
  delete cc;
  delete tf1; 
}

//
void GetPeakInfo(double par[2], vector<double> &x, vector<double> &y, map<double, double> &m_p)
{
  //par: calibration par
  //x: channel info
  //y: counts info
  //m_p: peak info from mca spectrum
  //histogram
  TCanvas *cc = new TCanvas("cc", "cc", 0, 0, 1200, 800);
  cc->Divide(2, 2);
  int bin_number = x.size();
  TH1D* h0 = new TH1D("h0", "raw data", bin_number, x[0], x[bin_number-1]);
  TH1D* h1 = new TH1D("h1", "calibration data", bin_number, par[0]+par[1]*x[0], par[0]+par[1]*x[bin_number-1]);
  for(int i=0;i<bin_number;i++){
    h0->SetBinContent(x[i], y[i]);
    h1->SetBinContent(x[i], y[i]);
  }
  cc->cd(1);
  h0->GetXaxis()->SetTitle("Channel");
  h0->GetYaxis()->SetTitle("Count");
  h0->Draw();
  cc->cd(2);
  h1->GetXaxis()->SetTitle("Energy[keV]");
  h1->GetYaxis()->SetTitle("Count");
  h1->Draw();

  //find peak
  cc->cd(3);
  TSpectrum *spec = new TSpectrum(20);
  int nfound = spec->Search(h0, 1, "", 0.01);
  TH1 *hb = spec->Background(h0, 20, "Compton BackSmoothing13 same");

  TF1 *tf_gaus = new TF1("tf_gaus", "gaus");

  double *xpeaks;
  xpeaks = spec->GetPositionX();
  for(int i=0;i<nfound;i++){
    double xp = xpeaks[i];
    if((par[0]+par[1]*xp) < 1.5)  continue;

    h0->Fit(tf_gaus, "", "", 0.97*xp, 1.03*xp);
    cout << "xp  " << xp << "  " << tf_gaus->GetParameter(1) << endl;

    map<double, double> m_single_peak;
    map<double, double>::iterator it;
    for(int j=0.9*xp;j<1.1*xp;j++){
      m_single_peak.insert(pair<double, double>(h0->GetBinCenter(j), h0->GetBinContent(j)));
    }
    for(it=m_single_peak.begin();it!=m_single_peak.end();++it){
      cout << it->first << " ==> " << it->second << endl;
    }
    map<double, double> m_find_max_x;
    for(int j=xp-2;j<=xp+2;j++){
      m_find_max_x.insert(pair<double, double>(h0->GetBinContent(j), h0->GetBinCenter(j)));
    }
    /*
    for(it=m_find_max_x.begin();it!=m_find_max_x.end();++it){
      cout << it->first << " ==> " << it->second << endl;
    }
    */
    map<double, double>::reverse_iterator rit = m_find_max_x.rbegin();
    cout << " !! " << rit->first << " " << rit->second << endl;

    double x_left_stop = 0;
    double x_right_stop = 0;
    //move left
    it = m_single_peak.find(rit->second);
    double y_max = it->second;
    cout << "y_max  " << y_max << endl;
    while(1){
      --it;
      if(it->second < y_max){
        y_max = it->second;
      }
      else{
        ++it;
        x_left_stop = it->first;
        break;
      }
    }
    //move right
    it = m_single_peak.find(rit->second);
    y_max = it->second;
    cout << "y_max  " << y_max << endl;
    while(1){
      ++it;
      if(it->second < y_max){
        y_max = it->second;
      }
      else{
        --it;
        x_right_stop = it->first;
        break;
      }
    }
    cout << "x_left_stop = " << x_left_stop << endl;
    cout << "x_right_stop = " << x_right_stop << endl;

    double area = h1->Integral(x_left_stop, x_right_stop) - hb->Integral(x_left_stop, x_right_stop);
    m_p.insert(pair<double, double>(par[0]+par[1]*xp, area));
  }
}

//
void CaliDetectorEff(double par1[10], double par2[10], double par3[10], map<double, double> &m_p, map<double, double> &m_pp)
{
  //par1, par2, par3: detector eff par
  //m_p: peak info from mca spectrum
  //m_pp: peak info from mca spectrum after detector efficiency calibration
  double eff = 0.;
  map<double, double>::iterator it;
  for(it=m_p.begin();it!=m_p.end();++it){
    eff = 0.;
    for(int i=0;i<10;i++){
      if(it->first<EFFENERGY1){
        eff += par1[i]*pow(it->first, (double)i);
      }
      else if(it->first>EFFENERGY2){
        eff += par3[i]*pow(it->first, (double)i);
      }else{
        eff += par2[i]*pow(it->first, (double)i);
      }
    }
    cout << "eff = " << eff << endl;
    m_pp.insert(pair<double, double>(it->first, it->second/eff));
  }
}

//
void CaliAttenEff(map<double, double> &m_p, map<double, double> &m_pp)
{
  //m_p: peak info from mca spectrum after detector efficiency calibration
  //m_pp: peak info from mca spectrum after detector&atten efficiency calibration
  ifstream fi;
  fi.open("gamma_attenuation_data.txt");
  if(!fi){
    cout << "cannot open gamma_attenuation_data.txt file !" << endl;
    return;
  }
  vector<double> v_x_ray_e;
  vector<double> v_x_ray_eff;

  double a, b;
  while(1){
    fi >> a >> b;
    if(!fi.good()) break;
    v_x_ray_e.push_back(a);
    v_x_ray_eff.push_back(b);
  }
  fi.close();

  int size = v_x_ray_e.size();
  TGraph *gr = new TGraph(size);
  for(int i=0;i<size;i++){
    gr->SetPoint(i, v_x_ray_e[i], v_x_ray_eff[i]);
  }

  map<double, double>::iterator it;
  for(it=m_p.begin();it!=m_p.end();++it){
    cout << "atten = " << gr->Eval(it->first) << endl;
    m_pp.insert(pair<double, double>(it->first, it->second/gr->Eval(it->first)));
  }
}

//
void GetElementShow(vector<string> &v_e, map<string, map<double, double>> &m_p, map<double, double> &m_pp, map<string, map<double, double>> &m_ppp, map<string, map<double, double>> &m_q)
{
  //v_e: the elements wants to caliculate
  //m_p: x-ray data from toi
  //m_pp: peak info from mca spectrum after detector&atten efficiency calibration
  //m_ppp: elements show in the mca spectrum, like <20Ca, <peak energy, peak area>, peak energy&area both from mca spectrum
  //m_q: associate peak energy from toi & peak from mca spectrum
  map<double, double> mm_p(m_pp);

  map<string, map<double, double>>::iterator itt1;
  map<double, double>::iterator itt2;

  int size = v_e.size();
  for(int i=0;i<size;i++){
    itt1 = m_p.find(v_e[i]);
    cout << "// analysis " << v_e[i] << " element." << endl;
    m_ppp.insert(pair<string, map<double, double>>(v_e[i], map<double, double>()));
    m_q.insert(pair<string, map<double, double>>(v_e[i], map<double, double>()));
    
    for(map<double, double>::iterator ittt=mm_p.begin();ittt!=mm_p.end();++ittt){
      cout << "//// " << ittt->first << endl;
    }

    for(itt2=itt1->second.begin();itt2!=itt1->second.end();++itt2){
      cout << "energy /// " << itt2->first << endl;
      map<double, double>::iterator ittt=mm_p.begin();
      double a = ittt->first;
      double d = abs(a-itt2->first);
      cout << "d = " << d << endl;
      while(1){
        ++ittt;
        a = ittt->first;
        if(d < abs(a-itt2->first)){
          --ittt;
          cout << "<< find peak " << ittt->first << endl;
          break;
        }
        else{
          d = abs(a-itt2->first);
          cout << "d = " << d << endl;
        }
      }
      //check real distance
      if(abs(ittt->first-itt2->first)/itt2->first > 0.1 || abs(d)>0.2){
        cout << "<< " << (ittt->first-itt2->first)/itt2->first << endl;
        cout << "!! can not find " << v_e[i] << " " << itt2->first << " keV X-Ray" << endl;
        continue;
      }
      //save and pop that peak
      m_ppp[v_e[i]].insert(pair<double, double>(ittt->first, ittt->second));
      m_q[v_e[i]].insert(pair<double, double>(ittt->first, itt2->first));
      mm_p.erase(ittt->first);
    }
  }

  //cout result
  for(map<string, map<double, double>>::iterator itttt=m_ppp.begin();itttt!=m_ppp.end();++itttt){
    cout << "<<<< element  " << itttt->first;
    for(itt2=itttt->second.begin();itt2!=itttt->second.end();++itt2){
      cout << "  " << itt2->first << "  " << itt2->second << endl;
    }
    cout << endl;
  }
  for(map<string, map<double, double>>::iterator itttt=m_q.begin();itttt!=m_q.end();++itttt){
    cout << "<<<< element  " << itttt->first;
    for(itt2=itttt->second.begin();itt2!=itttt->second.end();++itt2){
      cout << "  toi " << itt2->second << "keV &  mca " << itt2->first;
    }
    cout << endl;
  }

}

//
void CheckElementShow(map<string, map<double, double>> &m_p, map<string, map<double, double>> &m_pp, map<string, map<double, double>> &m_ppp, map<string, map<double, double>> &m_q)
{
  //m_p: x-ray data from toi
  //m_pp: elements show in the mca spectrum, like <20Ca, <peak energy, peak area>, peak energy&area both from mca spectrum
  //m_ppp: associate peak energy from toi & peak from mca spectrum, mca ==> toi
  //m_q: elements show in the mca spectrum after checked if more than one peak found, like <20Ca, <peak energy, peak area>, peak energy&area both from mca spectrum
  //peak energy from toi
  map<string, map<double, double>>::iterator it1;
  map<string, map<double, double>>::iterator it2;
  map<double, double>::iterator itt1;
  map<double, double>::iterator itt2;
  for(it1=m_pp.begin();it1!=m_pp.end();++it1){
    if(m_pp[it1->first].size() == 0){
      continue;
    }
    if(m_pp[it1->first].size() == 1){
      cout << ">>> only 1 peak found in " << it1->first << " go continue" << endl;
      itt1=m_ppp[it1->first].begin();
      itt2=m_pp[it1->first].begin();
      m_q[it1->first].insert(pair<double, double>(m_ppp[it1->first][itt1->first], itt2->second));
      continue;
    }
    else{//if more than 1 peak, choose the strongest one to make this easier.
         //if toi and mca show different strongest peak, than what ???
      for(itt1=m_p[it1->first].begin();itt1!=m_p[it1->first].end();++itt1){
        cout << ">>><< toi  " << it1->first << "  " << itt1->first << "  " << itt1->second << "  " << endl;
      }
      map<double, double> m_m;
      m_m.clear();
      for(itt1=m_pp[it1->first].begin();itt1!=m_pp[it1->first].end();++itt1){
        cout << ">>><< mca  " << it1->first << "  " << itt1->first << "  " << itt1->second << "  " << endl;
        m_m.insert(pair<double, double>(itt1->second, itt1->first));
      }
      map<double, double>::reverse_iterator rit = m_m.rbegin();
      cout << it1->first << " " << rit->second << " " << rit->first  << endl;
      m_q[it1->first].insert(pair<double, double>(m_ppp[it1->first][rit->second], rit->first));

      /*
      m_pp[it1->first].clear();
      m_pp[it1->first].insert(pair<double, double>(rit->second, rit->first));

      itt1=m_pp[it1->first].begin();
      cout << ">>>>> mca  " << it1->first << "  " << itt1->first << "  " << itt1->second << "  " << endl;
      */
    }
  }

  for(it2=m_q.begin();it2!=m_q.end();++it2){
    itt1=m_q[it2->first].begin();
    cout << ",,, " << it2->first << " " << itt1->first << " " << itt1->second << endl;
  }

  return;
}

//
void GetElementPercent(map<string, map<double, double>> &m_p, map<string, map<double, double>> &m_pp, map<string, double> &m_q)
{
  //m_p: x-ray data from toi
  //m_pp: elements show in the mca spectrum after checked if more than one peak found, like <20Ca, <peak energy, peak area>, peak energy&area both from mca spectrum 
  //m_q: element percent result

  map<string, map<double, double>>::iterator it;
  double sum = 0.;
  for(it=m_pp.begin();it!=m_pp.end();++it){
    map<double, double>::iterator itt = it->second.begin();
    sum += (itt->second)/(m_p[it->first][itt->first])*100;
  }
  cout << "sum = " << sum << endl;
  for(it=m_pp.begin();it!=m_pp.end();++it){
    map<double, double>::iterator itt = it->second.begin();
    m_q[it->first] = (itt->second)/(m_p[it->first][itt->first])*100/sum;
  }

  for(map<string, double>::iterator ittt=m_q.begin();ittt!=m_q.end();++ittt){
    cout << "!!!??? " << ittt->first << " " << ittt->second << endl;
  }
}