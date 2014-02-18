#include <stdio.h>

int regulate(int conc_count, double concs[conc_count],
		 int rows, int cols, double weights[rows][cols],
		 int tf_count, int tf_genes[tf_count],
		 int p_count, int p_genes[p_count],
		 int delta, double e_total){
  double updates[rows];
  double gene_count = rows - p_count;
  double signal;
  int below_zero = 0;
  int debug = 0;
  int idx;

  /*compute dot product of weight column and concs for update array*/
  for(int i = 0; i < cols; i++){
    signal = 0;
    for (int j = 0; j < rows; j++){
      signal += (concs[j] * weights[j][i]);
    }    
    signal = signal / gene_count;
    updates[i] = signal;
  }
  
  /*update the TF concentrations (scaled)*/
  for(int i = 0; i < tf_count; i++){
    idx = tf_genes[i];
    signal = (concs[idx] * updates[idx]) * delta;
    concs[idx] += signal;
    if(concs[idx] < 0){below_zero = 1;}
    if(concs[idx] < 1e-10){concs[idx] = 1e-10;}
  }

  /*update the P concentrations (direct)*/
  for(int i = 0; i < p_count; i++){
    idx = p_genes[i];
    concs[idx] += delta * updates[idx];
    if(concs[idx] < 1e-10){concs[idx] = 1e-10;}
  }

  /*calculate tf total scale for e_inputs and normalise tf genes*/
  double tf_total = 0;
  for(int i = 0; i < tf_count; i++){
    /*restart if TF conc drops below*/
    tf_total += concs[tf_genes[i]];
  }
  for(int i = 0; i < tf_count; i++){
    concs[tf_genes[i]] *= 1 - e_total;
    concs[tf_genes[i]] = concs[tf_genes[i]] / tf_total;
  } 

  /*calculate p total and normalise p genes*/
  double p_total = 0;
  for(int i = 0; i < p_count; i++){
    p_total += concs[p_genes[i]];
  }
  for(int i = 0; i < p_count; i++){
    concs[p_genes[i]] = concs[p_genes[i]] / p_total;
  } 

  return below_zero;
}

int regulateAsync(int conc_count, double concs[conc_count],
	     int rows, int cols, double weights[rows][cols],
	     int tf_count, int tf_genes[tf_count],
	     int p_count, int p_genes[p_count],
	     int delta, double e_total){
  double updates[rows];
  double gene_count = rows - p_count;
  double signal;
  int below_zero = 0;
  int debug = 0;

  /*compute dot product of weight column and concs for update array*/
  for(int i = 0; i < cols; i++){
    signal = 0;
    for (int j = 0; j < rows; j++){
      signal += (concs[j] * weights[j][i]);
    }    
    updates[i] = signal / gene_count;

    /*update the TF concentrations (scaled)*/
    for(int x = 0; x < tf_count; x++){
      if(i==tf_genes[x]){
	concs[tf_genes[x]] += (concs[tf_genes[x]] * updates[i]) * delta;
	/*restart if TF conc drops below*/
	if(concs[tf_genes[x]] < 0){below_zero = 1;}
	if(concs[tf_genes[x]] < 1e-10){concs[tf_genes[x]] = 1e-10;}
      }   
    }
    /*update the P concentrations (direct)*/
    for(int x = 0; x < p_count; x++){
      if(i==p_genes[x]){
	concs[p_genes[x]] += delta * updates[i];
	if(concs[p_genes[x]] < 1e-10){concs[p_genes[x]] = 1e-10;}
      }
    }
  }

  /*calculate tf total scale for e_inputs and normalise tf genes*/
  double tf_total = 0;
  for(int i = 0; i < tf_count; i++){
    tf_total += concs[tf_genes[i]];
  }
  for(int i = 0; i < tf_count; i++){
    concs[tf_genes[i]] *= 1 - e_total;
    concs[tf_genes[i]] = concs[tf_genes[i]] / tf_total;
  } 

  /*calculate p total and normalise p genes*/
  double p_total = 0;
  for(int i = 0; i < p_count; i++){
    p_total += concs[p_genes[i]];
  }

  for(int i = 0; i < p_count; i++){
    concs[p_genes[i]] = concs[p_genes[i]] / p_total;
  } 
  return below_zero;
}

