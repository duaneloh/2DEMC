/*
 *  EMC.c
 *  
 *  Created by N. Duane Loh 
 *  Last update: Tue Jun 10 14:01:37 SGT 2014
 *
 *	To compile :
 *		gcc -03 EMC.c -lm -o EMC
 *
 *	To use (1 argument) :
 *		./EMC num_of_iterations
 *	
 *	Needs: 
 *		data.dat (header format on line 94 of this file)
 *
 *	Output:
 *		image.dat
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

#define NEG_MIN (-1.E-308)
#define NUM_BACKGROUND (0)
#define LOW_PROBCUTOFF (1.E-308)
#define CHOICEUPDATEFREQ (100)
#define MINTOMOWT (1.E-308)
#define BGSUBITER (10)


struct global_var
{
	int len, tomo_len, conf_size, num_data, num_rlzt, num_background, present_conf, num_conf, num_tomo, iter, total_iter, choice_update;
	double err, mean_total_intens, mutual_info, m_info_given_choice;
	unsigned short **single_photon_pos, **multi_photon, *single_len, *multi_len;
	double *tomograms_in, *tomograms_out, *background_tomo, *tomograms_temp, *conf_in, *conf_out, *conf_out_wts, *p_rlzt, *cond_prob;
}G;

void setup();
void initialize_to_solution();
void initialize_to_random();
void free_mem();

void expand();
double maximize();
void compress();
void emc();
void print_recon();
void print_hitRates();

int main(int argc, char* argv[])
{
	//Clean up argc input.
	G.choice_update = 0;
	struct timeval tv1, tv2;
	G.num_conf = 2;
	G.num_background = 1;
	int i;

	if(argc == 1)
	{
		printf("Usage instructions::\n");
		printf("./EMC \n\t-n <num_iterations>");
		printf("\n\t-n 100 => iterate for 100 iterations.\n");
		return 0;
	}
	else if(argc > 1)
	{
		for(i = 0; i < argc; i++)
		{
			if(argv[i][0] == '-')
			{
				switch(argv[i][1])
				{
						case 'n':
							sscanf(argv[i+1], "%d", &G.total_iter);
							printf("number of iterations set to: %d\n", G.total_iter);
							break;
						/*
						You found the Easter egg! You can reconstruction multiple conformations if you used this flag!
						*/
						case 'c':
							sscanf(argv[i+1], "%d", &G.num_conf);
							G.num_conf += G.num_background;
							printf("num. of conformations set to: %d\n", G.num_conf);
							break;
				}
			}	
		}	
	}
	
	srand(time(0));
	setup();
	
	//Compute mutual information I(omega,K)|_W and write to logfile.
	initialize_to_solution();
	emc();
	FILE *fp;
	fp = fopen("EMC.log", "w");
	fprintf(fp, "num_data = %d \t mean_total_intens = %lf \t max_m_info = %lf \n", G.num_data, G.mean_total_intens, G.mutual_info);
	fclose(fp);
	fprintf(stderr, "max m_info = %lf \n", G.mutual_info);

	//Begin ab-initio reconstruction
	initialize_to_random();
    for(G.iter = 0; G.iter < G.total_iter; ++G.iter)
	{
		if(G.iter % CHOICEUPDATEFREQ == 0 && G.iter > 0)
			G.choice_update = 1;
		fp = fopen("EMC.log", "a");
		gettimeofday(&tv1, NULL);
		emc();
		gettimeofday(&tv2, NULL);
		//print_recon();
		fprintf(fp, "iteration = %d    error = %lf    m_info = %lf	 iter_time = %lf \n", G.iter + 1, G.err, G.mutual_info, (tv2.tv_sec-tv1.tv_sec) + (tv2.tv_usec-tv1.tv_usec)/(1000000.0F));
		fclose(fp);
		fprintf(stderr, "iteration = %d    error = %lf    m_info = %lf	 iter_time = %lf \n", G.iter + 1, G.err, G.mutual_info, (tv2.tv_sec-tv1.tv_sec) + (tv2.tv_usec-tv1.tv_usec)/(1000000.0F));
		G.choice_update = 0;
	}
    print_recon();
    print_hitRates();
	free_mem();
	return 0;
}

void setup()
{
	int r, c, t, d, i ;
	
	FILE *fptr;
	fptr = fopen("data.dat","rb");
	if(!fptr)
		fprintf(stderr, "Where is data.dat??\n"), exit(1);

	//fscanf(fptr, "%d  %d  %lf  %d  %d %lf ", &G.num_data, &G.present_conf, &G.mean_total_intens, &G.len, &G.tomo_len, &G.m_info_given_choice);
	fread(&G.num_data, sizeof(G.num_data), 1, fptr);
	fread(&G.present_conf, sizeof(G.present_conf), 1, fptr);
	fread(&G.mean_total_intens, sizeof(G.mean_total_intens), 1, fptr);
	fread(&G.len, sizeof(G.len), 1, fptr);
	fread(&G.tomo_len, sizeof(G.tomo_len), 1, fptr);
	fread(&G.m_info_given_choice, sizeof(G.m_info_given_choice), 1, fptr);
	
	G.conf_size= G.len * G.len ;
	G.num_tomo = 4 ;
	G.num_rlzt = ((G.num_conf-G.num_background) * (G.num_tomo) + G.num_background);
	fprintf(stderr, "Photon data has %d conformation(s).\nReconstructing with %d conformation(s) which includes %d background.\nHas %d data realizations.\n", G.present_conf, G.num_conf, G.num_background, G.num_rlzt);
	G.single_photon_pos = malloc(G.num_data * sizeof(*G.single_photon_pos));
	G.multi_photon = malloc(G.num_data * sizeof(*G.multi_photon));
	G.single_len = malloc(G.num_data * sizeof(*G.single_len));
	G.multi_len = malloc(G.num_data * sizeof(*G.multi_len));
	
	for(d = 0; d < G.num_data; ++d)
	{
		fread(&(G.single_len[d]), sizeof(*G.single_len),1,fptr);
//		fscanf(fptr, "%d ", &G.single_len[d]);
		G.single_photon_pos[d] = malloc(G.single_len[d] * sizeof(**G.single_photon_pos));
		fread((G.single_photon_pos[d]), sizeof(**G.single_photon_pos), G.single_len[d], fptr);
//		for(i = 0; i < G.single_len[d]; ++i)
//			fscanf(fptr, "%d ", &G.single_photon_pos[d][i]);

		fread(&(G.multi_len[d]), sizeof(*G.multi_len),1,fptr);
//		fscanf(fptr, "%d ", &G.multi_len[d]);
		G.multi_photon[d] = malloc(G.multi_len[d] * sizeof(**G.multi_photon));
		fread((G.multi_photon[d]), sizeof(**G.multi_photon), G.multi_len[d], fptr);	
//		for(i = 0; i < G.multi_len[d]; i+=2)
//			fscanf(fptr, "%d %d ", &G.multi_photon[d][i], &G.multi_photon[d][i+1]);
	}
	
	fclose(fptr);
	
	G.conf_in = malloc(G.num_conf * G.conf_size * sizeof(*G.conf_in));
	G.conf_out = malloc(G.num_conf * G.conf_size * sizeof(*G.conf_out));
	G.conf_out_wts = malloc(G.num_conf * G.conf_size * sizeof(*G.conf_out_wts));
	
	G.background_tomo = malloc(G.conf_size * sizeof(*G.background_tomo));
	G.tomograms_in = malloc(G.num_rlzt * G.conf_size * sizeof(*G.tomograms_in));
	G.tomograms_out = malloc(G.num_rlzt * G.conf_size * sizeof(*G.tomograms_out));
	G.tomograms_temp = malloc(G.num_rlzt * G.conf_size * sizeof(*G.tomograms_temp));
	G.p_rlzt = malloc(G.num_rlzt * sizeof(G.p_rlzt));
	G.cond_prob = malloc(G.num_data*G.num_rlzt * sizeof(*G.cond_prob));
	
	for(i = 0; i < G.num_rlzt; i++)
		G.p_rlzt[i] = 1./G.num_rlzt;
}

void initialize_to_solution() 
{
	int md, r, c, p;
	double scale_factor, total_contrast;
	FILE *fptr; 
	fptr = fopen("contrast.dat", "r");
	if(!fptr)
	{
		fprintf(stderr, "Cannot open contrast file!\n");
		exit(1); 
	}
	
	fscanf(fptr, "%d %d ", &G.len, &G.present_conf);
	if(G.conf_size != G.len * G.len)
		fprintf(stderr, "ImageSize in contrast.dat (%d) not equals that of data.dat (%d)!\n", G.len*G.len, G.conf_size), exit(1);
	
	for(md = 0; md < G.num_conf; ++md)
	{
		total_contrast = 0.;
		for(r = 0; r < G.len; ++r)
		for(c = 0; c < G.len; ++c)
		{
			p = (md * G.conf_size) + (r * G.len) + c;
			fscanf(fptr, "%lf ", &G.conf_in[p]);
			total_contrast += G.conf_in[p]; 
		}
		
		scale_factor = G.mean_total_intens / total_contrast;
		for(p = (md * G.conf_size); p < (md * G.conf_size) + G.conf_size; ++p)
			G.conf_in[p] *= scale_factor;
	}
	fclose(fptr);	
}

void initialize_to_random()
{
	int md, p, q;
    FILE *fp;
	double scale_factor, total_contrast, tempDoub;

	for(md = 0; md < G.num_conf; ++md)
	{
		total_contrast = 0.;
        if(md < G.num_conf-G.num_background)
        {
            for(p = (md * G.conf_size); p < (md * G.conf_size) + G.conf_size; ++p)
            {
                G.conf_in[p] = (((double) rand()) / RAND_MAX);
                total_contrast += G.conf_in[p];
            }
            
            scale_factor = G.mean_total_intens / total_contrast;
            for(p = (md * G.conf_size); p < (md * G.conf_size) + G.conf_size; ++p)
                G.conf_in[p] *= scale_factor;
        }
        else
        {
            fp = fopen("background.dat", "r");
        
            for(p = (md * G.conf_size); p < (md * G.conf_size) + G.conf_size; ++p)
            {
                fscanf(fp, "%lf ", &tempDoub);
                G.conf_in[p] = (((double) rand()) / RAND_MAX);
                //G.conf_in[p] = tempDoub;
                total_contrast += G.conf_in[p];
            }
            
            fclose(fp);
            scale_factor = G.mean_total_intens / total_contrast;
            for(p = (md * G.conf_size); p < (md * G.conf_size) + G.conf_size; ++p)
               G.conf_in[p] *= scale_factor;
        }
        
	}
    
    
}

void emc()
{
	int i;
	double temp_double;
	expand();
	G.mutual_info = maximize();
	G.mutual_info = 1. - (G.mutual_info / G.m_info_given_choice);
	compress();
	
	G.err = 0.; 
	for(i = 0; i < G.num_conf * G.conf_size; i++)
		{
		temp_double = G.conf_in[i] - G.conf_out[i];
		G.err += temp_double * temp_double;
		G.conf_in[i] = G.conf_out[i];
		}
	G.err = sqrt(G.err / (G.num_conf * G.conf_size)) / (G.mean_total_intens / (G.num_conf * G.conf_size));
}

void expand()
{
	int r, c, i, j, k, l, m, md;
	
	for(md = 0; md < G.num_conf; md++)
	{
		if(md < G.num_conf-G.num_background)
		{
			l = (md * G.num_tomo * G.conf_size);
			m = (md * G.conf_size);

			for(i = 0, r = 0; r < G.len; r++)
			for(c = 0; c < G.len; c++, i++)
				G.tomograms_in[l + i] = G.conf_in[m + (r * G.len) + c];

			for(i = 0, c = 0; c < G.len; c++)
			for(r = G.len - 1; r >= 0; r--, i++)
				G.tomograms_in[l + (G.conf_size) + i] = G.conf_in[m + (r * G.len) + c];

			for(i = 0, r = G.len - 1; r >= 0; r--)
			for(c = G.len - 1; c >= 0; c--, i++)
				G.tomograms_in[l + (2 * G.conf_size) + i] = G.conf_in[m + (r * G.len) + c];

			for(i = 0, c = G.len - 1; c >= 0; c--)
			for(r = 0; r < G.len; r++, i++)
				G.tomograms_in[l + (3 * G.conf_size) + i] = G.conf_in[m + (r * G.len) + c];
/*			
			for(i = 0, r = 0; r < G.len; r++)
			for(c = G.len - 1; c >= 0; c--, i++)
				G.tomograms_in[l + (4 * G.conf_size) + i] = G.conf_in[m + (r * G.len) + c];
			
			for(i = 0, c = 0; c < G.len; c++)
			for(r = 0; r < G.len; r++, i++)
				G.tomograms_in[l + (5 * G.conf_size) + i] = G.conf_in[m + (r * G.len) + c];
			
			for(i = 0, r = G.len - 1; r >= 0; r--)
			for(c = 0; c < G.len; c++, i++)
				G.tomograms_in[l + (6 * G.conf_size) + i] = G.conf_in[m + (r * G.len) + c];			
			
			for(i = 0, c = G.len - 1; c >= 0; c--)
			for(r = G.len - 1; r >= 0; r--, i++)
				G.tomograms_in[l + (7 * G.conf_size) + i] = G.conf_in[m + (r * G.len) + c];
*/			
		}
		else
		{
			l = (((G.num_conf-G.num_background) * G.num_tomo) + (md-(G.num_conf-G.num_background))) * G.conf_size;
			m = (md * G.conf_size);
			for(i = 0, r = 0; r < G.len; r++)
			for(c = 0; c < G.len; c++, i++)
				G.tomograms_in[l + i] = G.conf_in[m + (r * G.len) + c];
		}
	}
	
}

double maximize()
{
	int d, t, md, i, j, k;
	//todo: normalize probability of background.
	double p_tot, p_tomo_given_data, max_exp;
	double info = 0.;
	double s[G.num_rlzt], p[G.num_rlzt], p_choice_temp[G.num_rlzt], tomo_temp1[G.num_rlzt];
	
	for(md = 0; md < G.num_conf; md++)
	{
		if(md < G.num_conf-G.num_background)
		{
			for(t = 0; t < G.num_tomo; t++)
			{
				k = (md * G.num_tomo) + t;
				s[k] = 0.;
				p_choice_temp[k] = 0.;
				tomo_temp1[k] = log(G.p_rlzt[k]);
				for(i = 0; i < G.conf_size; i++)
				{
					j = (md * G.num_tomo * G.conf_size) + (t * G.conf_size) + i;
					G.tomograms_out[j]  = 0.;
					tomo_temp1[k] -= G.tomograms_in[j];
					G.tomograms_temp[j] = (G.tomograms_in[j] > 0) ? log(G.tomograms_in[j]) : NEG_MIN;
				}
			}
		}
		else
		{
			k = (G.num_conf-G.num_background)*G.num_tomo + (md-(G.num_conf-G.num_background));
			s[k] = 0.;
			p_choice_temp[k] = 0.;
			tomo_temp1[k] = log(G.p_rlzt[k]);
			t = G.conf_size * k;
			for(i = 0; i < G.conf_size; i++)
			{
				j = t + i;
				G.tomograms_out[j]  = 0.;
				tomo_temp1[k] -= G.tomograms_in[j];
				G.tomograms_temp[j] = (G.tomograms_in[j] > 0) ? log(G.tomograms_in[j]) : NEG_MIN;
			}
		}
	}
			
	for(d = 0; d < G.num_data; d++)
	{
		p_tot = 0.;
		max_exp = -100. * G.conf_size;
		//TODO: only operate on positive photon counts. Negative ones get incremented but not considered in conditional probability
		for(md = 0; md < G.num_conf; md++)
		{
			if(md < G.num_conf-G.num_background)
			{
				for(t = 0; t < G.num_tomo; t++)
				{
					k = (md * G.num_tomo) + t;
					p[k] = tomo_temp1[k];
					j = G.conf_size * k;
					for(i = 0; i < G.single_len[d]; ++i)
						p[k] += G.tomograms_temp[j + G.single_photon_pos[d][i]];
					
					for(i = 0; i < G.multi_len[d]; i += 2)
						p[k] += G.multi_photon[d][i+1] * G.tomograms_temp[j + G.multi_photon[d][i]];
					
					if(p[k] > max_exp)
						max_exp = p[k];
				}
			}
			else
			{
				k = ((G.num_conf-G.num_background) * G.num_tomo) + (md - (G.num_conf-G.num_background));
				p[k] = tomo_temp1[k];
				j = k * G.conf_size;
				for(i = 0; i < G.single_len[d]; ++i)
					p[k] += G.tomograms_temp[j + G.single_photon_pos[d][i]];
				
				for(i = 0; i < G.multi_len[d]; i += 2)
					p[k] += G.multi_photon[d][i+1] * G.tomograms_temp[j + G.multi_photon[d][i]];
					
				if(p[k] > max_exp)
					max_exp = p[k];
			}
		}
		
		for(md = 0; md < G.num_rlzt; md++)
		{
			p[md] = exp(p[md]-max_exp);
			p_tot += p[md];
		}

		for(md = 0; md < G.num_conf; md++)
		{
			if(md < G.num_conf-G.num_background)
			{
				for(t = 0; t < G.num_tomo; t++)
				{
					k = (md * G.num_tomo) + t;
					p[k] /= p_tot;
					G.cond_prob[d*G.num_rlzt + k] = p[k];
					s[k] += p[k];
					p_choice_temp[k] += p[k];
					if(p[k] > LOW_PROBCUTOFF) //TODO: Set a higher cutoff.
						info += p[k] * log(p[k] / G.p_rlzt[k]);
					else
						continue;

					j = k * G.conf_size;
					for(i = 0; i < G.single_len[d]; i++)
						G.tomograms_out[j + G.single_photon_pos[d][i]]  += p[k];
					for(i = 0; i < G.multi_len[d]; i += 2)
						G.tomograms_out[j + G.multi_photon[d][i]] += p[k] * G.multi_photon[d][i+1];
				}
			}
			else
			{
				k = ((G.num_conf-G.num_background) * G.num_tomo) + (md - (G.num_conf-G.num_background));
				p[k] /= p_tot;
				G.cond_prob[d*G.num_rlzt + k] = p[k];
				s[k] += p[k];
				p_choice_temp[k] += p[k];
				if(p[k] > LOW_PROBCUTOFF) //TODO: Set a higher cutoff.
					info += p[k] * log(p[k] / G.p_rlzt[k]);
				else
					continue;

				j = k * G.conf_size;				
				for(i = 0; i < G.single_len[d]; i++)
					G.tomograms_out[j + G.single_photon_pos[d][i] ]  += p[k];
				for(i = 0; i < G.multi_len[d]; i += 2)
					G.tomograms_out[j + G.multi_photon[d][i] ] += p[k] * G.multi_photon[d][i+1];
				
			}
		}
		
	}
	
	double tot_rlzt = 0.;
	for(md = 0; md < G.num_rlzt; md++)
	{
		if(s[md] < 1.E-308)
			continue;
		tot_rlzt += p_choice_temp[md];

		k = md * G.conf_size;
		for(i = 0; i < G.len; ++i)
		for(j = 0; j < G.len; ++j)
			G.tomograms_out[k + i*G.len + j] /= s[md];
	}
	if(G.choice_update == 1)
	{
		for(md = 0; md < G.num_rlzt; md++)
		{
			G.p_rlzt[md] = p_choice_temp[md]/tot_rlzt; 
			fprintf(stderr, "%lf ", G.p_rlzt[md]);
		}
		fprintf(stderr, "\n");
	}	
	
	info /= 1.*G.num_data;

	
/*
	//Bypass maximization step: used to test expansion and compression. 
	for(md = 0; md < G.num_rlzt; md++)
	{
		k = md * G.conf_size;
		for(i = 0; i < G.len; ++i)
		for(j = 0; j < G.len; ++j)
			G.tomograms_out[k + i*G.len + j] = G.tomograms_in[k + i*G.len + j];
	}
*/
	return info;
}

void compress()
{
	int r, c, i, j, k, l, m, md;

	for(i = 0; i < G.num_rlzt*G.conf_size; i++)
		{
		G.conf_out[i] = 0.;
		G.conf_out_wts[i] = 0.;
		}
	for(md = 0; md < G.num_conf; md++)
	{
		if(md < G.num_conf-G.num_background)
		{
			l = (md * G.num_tomo * G.conf_size);
			m = (md * G.conf_size);
			
			for(i = 0, r = 0; r < G.len; r++)
			for(c = 0; c < G.len; c++, i++)
			{
				if(G.tomograms_out[l+i] >= MINTOMOWT)
				{
					G.conf_out[m + (r*G.len) + c] += G.tomograms_out[l+i];
					G.conf_out_wts[m + (r*G.len) + c] += 1.;
				}
			}				
			
			for(i = 0, c = 0; c < G.len; c++)
			for(r = G.len - 1; r >= 0; r--, i++)
			{
				if(G.tomograms_out[l + (G.conf_size) + i] >= MINTOMOWT)
				{
					G.conf_out[m + (r * G.len) + c] += G.tomograms_out[l + (G.conf_size) + i];
					G.conf_out_wts[m + (r*G.len) + c] += 1.;
				}
			}
			
			for(i = 0, r = G.len - 1; r >= 0; r--)
			for(c = G.len - 1; c >= 0; c--, i++)
			{
				if(G.tomograms_out[l + (2*G.conf_size) + i] >= MINTOMOWT)
				{
					G.conf_out[m + (r * G.len) + c] += G.tomograms_out[l + (2*G.conf_size) + i];
					G.conf_out_wts[m + (r*G.len) + c] += 1.;
				}
			}
			
			for(i = 0, c = G.len - 1; c >= 0; c--)
			for(r = 0; r < G.len; r++, i++)
			{
				if(G.tomograms_out[l + (3*G.conf_size) + i] >= MINTOMOWT)
				{
					G.conf_out[m + (r * G.len) + c] += G.tomograms_out[l + (3*G.conf_size) + i];	
					G.conf_out_wts[m + (r*G.len) + c] += 1.;
				}
			}
/*
			for(i = 0, r = 0; r < G.len; r++)
			for(c = G.len - 1; c >= 0; c--, i++)
			{
				if(G.tomograms_out[l + (4*G.conf_size) + i] >= MINTOMOWT)
				{
					G.conf_out[m + (r * G.len) + c] += G.tomograms_out[l + (4*G.conf_size) + i];
					G.conf_out_wts[m + (r*G.len) + c] += 1.;
				}
			}
			
			for(i = 0, c = 0; c < G.len; c++)
			for(r = 0; r < G.len; r++, i++)
			{
				G.conf_out[m + (r * G.len) + c] += G.tomograms_out[l + (5*G.conf_size) + i];
				if(G.tomograms_out[l + (5*G.conf_size) + i] >= MINTOMOWT)
					G.conf_out_wts[m + (r*G.len) + c] += 1.;
			}
			
			for(i = 0, r = G.len - 1; r >= 0; r--)
			for(c = 0; c < G.len; c++, i++)
			{
				if(G.tomograms_out[l + (6*G.conf_size) + i] >= MINTOMOWT)
				{
					G.conf_out[m + (r * G.len) + c] += G.tomograms_out[l + (6*G.conf_size) + i];
					G.conf_out_wts[m + (r*G.len) + c] += 1.;
				}
			}
			
			for(i = 0, c = G.len - 1; c >= 0; c--)
			for(r = G.len - 1; r >= 0; r--, i++)
			{
				if(G.tomograms_out[l + (7*G.conf_size) + i] >= MINTOMOWT)
				{
					G.conf_out[m + (r * G.len) + c] += G.tomograms_out[l + (7*G.conf_size) + i];
					G.conf_out_wts[m + (r*G.len) + c] += 1.;
				}
			}
*/
		}
		else
		{
			l = (((G.num_conf-G.num_background) * G.num_tomo) + (md-(G.num_conf-G.num_background))) * G.conf_size;
			m = (md * G.conf_size);
	
			for(i = 0, r = 0; r < G.len; r++)
			for(c = 0; c < G.len; c++, i++)
				G.conf_out[m + (r * G.len) + c] += G.tomograms_out[l + i];
		}
	}

	double norm;	
	for(md = 0; md < G.num_conf; md++)
	{
		if(md < G.num_conf-G.num_background)
		{
			m = (md * G.conf_size);
			norm = 1.0 / G.num_tomo; //Assumes uniform orientation distribution.
			for(r = 0; r < G.len; r++)
			for(c = 0; c < G.len; c++)
				{
				if(G.conf_out_wts[m + (r*G.len) + c] > 1.E-7)
					G.conf_out[m + (r*G.len) + c] /= G.conf_out_wts[m + (r*G.len) + c];
				}
		}
	}

}	

void print_recon()
{
	int r, c, md, i, j, k;
	FILE *fp, *fp2;
	char output_name[100];
	sprintf(output_name, "recon%03d.dat", G.iter);
	fp = fopen(output_name, "w");
	for(i = 0, md = 0; md < G.num_conf; md++)
	{
		for(r = 0; r < G.len; r++)
		{	
			for(c = 0; c < G.len; c++, i++)
				fprintf(fp, "%lf ", G.conf_out[i]);
			
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	sprintf(output_name, "cond_prob%03d.dat", G.iter);
	fp = fopen(output_name, "w");
	for(i=0, k=0; i<G.num_data; i++)
	{
		for(j=0; j<G.num_rlzt; j++, k++)
			fprintf(fp, "%e ", G.cond_prob[k]);
		fprintf(fp, "\n");
	}
	fclose(fp);

}

void print_hitRates()
{
	int r, c, md, i, j, k;
	FILE *fp, *fp2;
	char input_name[100];
    char output_name[100];
    
	sprintf(input_name, "hidden_variables.dat");
	fp2 = fopen(input_name, "r");
    double ftemp, ftemp2, tot, hits, falseHits, blanks, falseBlanks;
    tot = hits = falseHits = blanks = falseBlanks = 0.;
    double tempPh, hitCount, hitPh, hitPhSq, blankCount, blankPh, blankPhSq;
    hitCount = hitPh = hitPhSq = blankCount = blankPh = blankPhSq = 0.;
    
	for(i=0, k=0; i<G.num_data; i++, k+=G.num_rlzt)
	{
        fscanf(fp2, "%d %lf %lf ", &r, &ftemp, &ftemp2);
        tempPh = G.single_len[i];
        for(j = 0; j < G.multi_len[i]; j += 2)
            tempPh += G.multi_photon[i][j+1] ;

        if(r != -1)
        {
            hitCount += 1.0;
            hitPh += tempPh;
            hitPhSq += tempPh*tempPh ;
        }
        else
        {
            blankCount += 1.0;
            blankPh += tempPh;
            blankPhSq += tempPh*tempPh ;
        }
        
        for(j=0; j < G.num_rlzt; j++)
        {
            tot += G.cond_prob[k+j];
            if(j<(G.num_conf-G.num_background)*G.num_tomo)
            {
                if(r != -1)
                    hits += G.cond_prob[k+j];
                else
                    falseHits += G.cond_prob[k+j];
            }
            else
            {
                if(r == -1)
                    blanks += G.cond_prob[k+j];
                else
                    falseBlanks += G.cond_prob[k+j];
            }
        }
	}
	fclose(fp2);
    
    if(hitCount > 0.)
    {
        hitPh /= hitCount;
        hitPhSq /= hitCount;
        hitPhSq = sqrt(hitPhSq - hitPh*hitPh);
    }

    if(blankCount > 0.)
    {
        blankPh /= blankCount;
        blankPhSq /= blankCount;
        blankPhSq = sqrt(blankPhSq - blankPh*blankPh);
    }
    
    
	sprintf(output_name, "hitFractions.dat");
	fp = fopen(output_name, "a");
    fprintf(fp, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", hitPh, hitPhSq, blankPh, blankPhSq, hits/tot, falseHits/tot, blanks/tot, falseBlanks/tot);
	fclose(fp);
    
}

void free_mem()
{
	int d;
	
	for(d = 0; d < G.num_data; ++d) 
		free(G.single_photon_pos[d]), free(G.multi_photon[d]);
	free(G.single_photon_pos), free(G.multi_photon);
	free(G.single_len), free(G.multi_len);
	
	free(G.conf_in); 
	free(G.conf_out);
	free(G.conf_out_wts);
	
	free(G.background_tomo);
	free(G.tomograms_in);
	free(G.tomograms_out);
	free(G.tomograms_temp);
	free(G.p_rlzt);
	free(G.cond_prob);
}
