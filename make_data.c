/*
 *  make_data.c
 *  
 *  Created by Duane Loh 
 *  Last update: Tue Jun 10 14:01:37 SGT 2014
 *
 *	To compile :
 *		gcc -03 make_data.c -lm -o make_data
 *
 *	Usage instructions:
 *		./make_data
 *	
 *	Needs:
 *		contrast.dat (header format on line 119 of this file)
 *
 *	Output:
 *		data.dat (header format on line 63 of this file)
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define PI (3.1415926535897932)
#define CALIB_CYCLE (1000)
#define K_DIV (1000)
#define PKW_INNER_CYCLE (1000)
#define DEFAULT_SCALE_DIV (200000.)
#define SCALE_MIN (1.0)
#define SCALE_MAX (1.0)
#define DEFAULT_HIT_RATE (1.)
#define FLUENCE_VAR (0.05)
#define TOMOFRAC (1.)

void initialize();
void expand_image();
void free_mem();
int poisson_dist(double);
double gaussian_dist(double, double);
double exp_dist(double); 
double fac(int);

int len, img_size, tomo_len, num_tomo, num_data, num_imgs, data_counter;
double half_len, mean_total_photons, mean_intens, mean_total_noise_var, avg_noise_counter;
double *img , *tomo, *noise, *avg_noise;
int *single_photons, *multi_photons;

int main(int argc, char* argv[])
{
	
	int rand_choice, im, t, r, c, pos, i, j, single_len, multi_len, curr_photon1, curr_photon2, curr_intens_pos, int_conformation_div, max_photon = 0;
	double *k_hist, *w_hist, *kw_hist;
	double photon_scaling, tot_signal_photons, rand_float, m_info_given_choice, m_info_given_choice_counter, hit_rate, conformation_div;
	
	num_data = 100;
	mean_total_photons = 100.;
	hit_rate = DEFAULT_HIT_RATE;
	mean_total_noise_var = 0.;
	conformation_div = DEFAULT_SCALE_DIV;
	
	if(argc == 1)
	{
		printf("Usage instructions::\n");
		printf("./make_data \n\t-s <%.1lf> -d <%d> -n <%.1lf> -h <%.1lf> -c <%d>; default values in angle brackets;\n", mean_total_photons, num_data, mean_total_noise_var, hit_rate, (int) conformation_div);
		printf("\n\t-s 500. => input signal photons per data set to 500.0\n");
		printf("\t-d 1000 => num of data set to 1000\n");
		printf("\t-n 100. => input noise level per data set to 100.0\n");
		printf("\t-h 0.4  => hit rate set to 40 percent\n");
		printf("\t-c 10	  => number of conformations set to 10\n");

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
						case 's':
							sscanf(argv[i+1], "%lf", &mean_total_photons);
							printf("signal photon level set to: %lf\n", mean_total_photons);
							break;
						case 'd':
							sscanf(argv[i+1], "%d", &num_data);
							printf("num data set to: %d\n", num_data);
							break;
						case 'n':
							sscanf(argv[i+1], "%lf", &mean_total_noise_var);
							printf("noise photon level set to: %lf\n", mean_total_noise_var);
							break;
						case 'h':
							sscanf(argv[i+1], "%lf", &hit_rate);
							printf("hit rate set to: %lf\n", hit_rate);
							break;
						case 'c':
							sscanf(argv[i+1], "%d", &int_conformation_div);
							conformation_div = (double) int_conformation_div;
							printf("num. of conformations set to: %d\n", (int) conformation_div);
							break;
				}
			}	
		}	
	}
	
	num_tomo = 4;
	srand(time(0));
	initialize();
	expand_image();

	FILE *fp, *fps;
	fp = fopen("data.dat" ,"w");
	fps = fopen("hidden_variables.dat", "w");
	
	double calib_signal_counter = 0.;	
	avg_noise_counter = 0.;
	tot_signal_photons = 0.;
	photon_scaling = 1.;
	for(data_counter = 0; data_counter < num_data+CALIB_CYCLE; ++data_counter)
	{
		single_len = 0; multi_len = 0;
		
		if(data_counter == CALIB_CYCLE)
		{
			calib_signal_counter /= ((double) tomo_len*tomo_len);
			photon_scaling = mean_total_photons/(tot_signal_photons/calib_signal_counter);
			printf("mean_tot_photons: %lf, \tphoton scaling: %lf\n", mean_total_photons, photon_scaling);
		}
		rand_choice = (((double)rand())/RAND_MAX < hit_rate)? ((int) (((double)rand()*num_tomo*num_imgs)/RAND_MAX)) : -1; 
		double fluence = gaussian_dist(photon_scaling, FLUENCE_VAR);
		if(fluence < 0.) fluence = 0.;

		//Mutual information computed only once when data_counter == CALIB_CYCLE. 
		//TODO: mutual information computation hasn't included model scaling.
		double total_prob = 0.;
		if(data_counter == CALIB_CYCLE)
		{
			max_photon *= 10;
			k_hist = malloc((max_photon+1) * sizeof(*k_hist));
			for(i = 0 ; i <= max_photon ; i++)	
				{k_hist[i] = 0.;}
			w_hist = malloc((K_DIV*(max_photon+1)) * sizeof(*w_hist));
			kw_hist = malloc((max_photon+1) * (K_DIV*(max_photon+1)) * sizeof(*kw_hist)); 	

			m_info_given_choice = 0., m_info_given_choice_counter = 0.;

			double k_hist_tot = 0., kw_hist_tot;
			for(i = 0; i <= max_photon; i++)	
				{k_hist_tot += k_hist[i];}
			for(i = 0; i <= max_photon; i++)	
				{k_hist[i] /= k_hist_tot;}

			kw_hist_tot = 0.;
			for(i = 0; i < (max_photon+1) * (K_DIV*(max_photon+1)); i++)
				{kw_hist[i] = 0.;}

			//TODO: change to random interp model?
			for(im = 0; im < num_imgs; im++)
			for(t = 0; t < img_size; t++)
			{
				pos = (im * num_tomo * img_size) + t;

				for(i = 0; i < PKW_INNER_CYCLE; i++)
				{
					curr_intens_pos = (int) (K_DIV * tomo[pos]);
					curr_photon1 = poisson_dist(fluence*tomo[pos]);
					curr_photon2 = poisson_dist(noise[r*len+c]);
					//curr_photon1 = (curr_photon1 > 0) ? curr_photon1 : 0;
					//curr_photon2 = (curr_photon2 > 0) ? curr_photon2 : 0;
					curr_photon1 += curr_photon2;
					if(curr_photon1 <= max_photon && curr_intens_pos <= K_DIV*(max_photon))
					{
						kw_hist[curr_photon1*(K_DIV*(max_photon+1)) + curr_intens_pos] += 1.;
						kw_hist_tot += 1.;	
					}
				}			
			}
				
			for(i = 0; i < (max_photon+1) * (K_DIV*(max_photon+1)); i++)
				{kw_hist[i] /= kw_hist_tot;}
			
			double temp;
			for(i = 0; i <= max_photon; i++)
			{
				temp = 0.;
				for(j = 0; j <= K_DIV*max_photon; j++)
					temp += kw_hist[i*K_DIV*(max_photon+1) + j];
				k_hist[i] = temp;
			}
			for(i = 0; i <= K_DIV*max_photon; i++)
			{
				temp = 0.;
				for(j = 0; j <= max_photon; j++)
					temp += kw_hist[j*(K_DIV*(max_photon+1)) + i];
				w_hist[i] = temp;
			}
			for(i = 0; i <= max_photon; i++)
			for(j = 0; j <= K_DIV*max_photon; j++)
			{
				temp = kw_hist[i*K_DIV*(max_photon+1) + j];
				if(temp > 0.)
					m_info_given_choice += temp * log(temp / (k_hist[i] * w_hist[j]));
				total_prob += temp;	
			}
					
			m_info_given_choice *= img_size;
			fprintf(fp, "%d %d %e %d %d %e\n", num_data, (int)(num_imgs*conformation_div), mean_total_photons, len, tomo_len, m_info_given_choice);
		}
				
		int r0, c0, r1, c1, tomo_offset;
		double v0, v1, v2, v3, v_interp, dc, dr, fc, fr;
		double curr_scale = floor( ((conformation_div*rand())/RAND_MAX) ) / conformation_div;
		double gamma = SCALE_MIN + (SCALE_MAX-SCALE_MIN)*curr_scale;
		if(gamma > SCALE_MAX)		{gamma = SCALE_MAX;}
		else if(gamma < SCALE_MIN)  {gamma = SCALE_MIN;}
		for(r = 0; r < tomo_len; r++)
		for(c = 0; c < tomo_len; c++)
		{
			t = r*len + c;
			if(rand_choice >= 0)
			{
				fr = gamma*(r-half_len) + half_len;
				fc = gamma*(c-half_len) + half_len;
				r0 = (int) fr; r1 = r0+1; dr = fr-r0;
				c0 = (int) fc; c1 = c0+1; dc = fc-c0;
				tomo_offset = (rand_choice * img_size);
				v0 = tomo[tomo_offset + r0*len + c0];
				v1 = tomo[tomo_offset + r0*len + c1];
				v2 = tomo[tomo_offset + r1*len + c0];
				v3 = tomo[tomo_offset + r1*len + c1];
				v_interp = fluence*(v0 + dc*dr*(v0-v1-v2+v3) + dc*(v1-v0) +dr*(v2-v0));
				
				curr_photon1 = poisson_dist(v_interp);
				//curr_photon2 = 0;
				curr_photon2 = poisson_dist(noise[r*len+c]);
				//curr_photon1 = (curr_photon1 > 0) ? curr_photon1 : 0;
				//curr_photon2 = (curr_photon2 > 0) ? curr_photon2 : 0;
			}
			else if(rand_choice == -1)
			{
				curr_photon1 = 0;
				curr_photon2 = poisson_dist(noise[r*len+c]);
				//curr_photon2 = (curr_photon2 > 0) ? curr_photon2 : 0;
				if(data_counter >= CALIB_CYCLE)
					avg_noise[r*len+c] += curr_photon2;
			}

			if(data_counter < CALIB_CYCLE)
			{
				if(rand_choice >= 0)
				{
					calib_signal_counter += 1.0;
					tot_signal_photons += curr_photon1;
				}
				
				if(curr_photon1+curr_photon2 > max_photon)
					max_photon = curr_photon1+curr_photon2;
					
				continue ;
			}

			curr_photon1 += curr_photon2;
											
			if(curr_photon1 == 1) 
				single_photons[single_len++] = t;
			else if(curr_photon1 > 1)
			{
				multi_photons[multi_len++] = t;
				multi_photons[multi_len++] = curr_photon1;
			}
		}
			
		if(data_counter < CALIB_CYCLE)
			continue ;

		if(rand_choice < 0)
			avg_noise_counter += 1.;
							
		fprintf(fp, "%d\n", single_len) ;
		for(t = 0 ; t < single_len ; ++t)
			fprintf(fp, "%d ", single_photons[t]) ;
		fprintf(fp, "\n") ;

		fprintf(fp, "%d\n", multi_len) ;
		for(t = 0 ; t < multi_len ; t+=2)
			fprintf(fp, "%d %d ", multi_photons[t], multi_photons[t+1]) ;
		fprintf(fp, "\n\n") ;
			
		fprintf(fps, "%d %lf %lf\n", rand_choice, gamma, fluence);
	}
	
	fprintf(fp, "\n");
	fclose(fp);
	fclose(fps);
		
	fp = fopen("background.dat", "w");
	for(r = 0; r < len; r++)
	{
		for(c = 0; c < len; c++)
		{
			fprintf(fp, "%lf ", avg_noise[r*len+c]/avg_noise_counter);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	
	free(k_hist), free(kw_hist);
	free_mem();
	return 0;
	
}
	
void initialize() 
{
	int im, r, c, p;
	double scale_factor, total_contrast;

	FILE *fptr; 
	fptr = fopen("contrast.dat", "r");
	if(!fptr)
	{
		fprintf(stderr, "Cannot open contrast file!\n");
		exit(1); 
	}
	
	fscanf(fptr, "%d %d ", &len, &num_imgs);
	tomo_len = TOMOFRAC * len;
	img_size= len * len;
	half_len = ( (len & 1) == 0 ) ? (len/2. - 0.5) : ( (len - 1)/2. );
	img = malloc( num_imgs * img_size * sizeof(*img) );
	noise = malloc( img_size * sizeof(*noise) );
	avg_noise = malloc(img_size * sizeof(*avg_noise));
	single_photons = malloc(img_size* sizeof(*single_photons));
	multi_photons = malloc(2 * img_size* sizeof(*multi_photons));
	mean_total_noise_var /= (double)img_size;
	
	for(im = 0; im < num_imgs; ++im)
	{
		total_contrast = 0.;
		for( r = 0; r < len; ++r)
		for( c = 0; c < len; ++c)
		{
			p = (im*img_size) + (r*len) + c;
			fscanf(fptr, "%lf ", &img[p]);
			total_contrast += img[p]; 
		}
		
		scale_factor = mean_total_photons / total_contrast;
		for(p = (im*img_size); p < (im*img_size) + img_size; ++p)
			img[p] *= scale_factor;
	}
	fclose(fptr) ;
	
	tomo = malloc(num_imgs * num_tomo * img_size * sizeof(*tomo));
	
	double noise_norm = 1., total_noise = 0., total_noise1 = 0.;	
	total_noise = 0.;
	for(r=0; r<len; r++)
	for(c=0; c<len; c++)
	{
		noise[r*len + c] = exp_dist(mean_total_noise_var);
		total_noise += poisson_dist(noise[r*len + c]);
		avg_noise[r*len + c] = 0.;
	}
	total_noise /= (double) img_size;
	noise_norm = mean_total_noise_var/total_noise;
}
	
	
void expand_image()
{
	int r, c, p, t, i, im ;
	for(im = 0; im < num_imgs; ++im)
	{
		//Rot0, Flip0
		for(i = 0, r = 0; r < len; ++r)
		for(c = 0; c < len; ++c, ++i)
		{
			p = (im * img_size) + (r * len) + c;
			t = (im * num_tomo * img_size) + i;
			tomo[t] = img[p];
		}
		
		//Rot90, Flip0
		for(i = 0, c = 0; c < len; ++c)
		for(r = len - 1; r >= 0; --r, ++i)
		{
			p = (im * img_size) + (r * len) + c;
			t = (im * num_tomo * img_size) + (img_size) + i;
			tomo[t] = img[p];
		}

		//Rot180, Flip0
		for(i = 0, r = len - 1; r >= 0; --r)
		for(c = len - 1; c >= 0; --c, ++i)
		{
			p = (im * img_size) + (r * len) + c;
			t = (im * num_tomo * img_size) + (2 * img_size) + i;
			tomo[t] = img[p];
		}
		
		//Rot270, Flip0
		for(i = 0, c = len - 1 ; c >= 0 ; --c)
		for(r = 0 ; r < len ; ++r, ++i)
		{
			p = (im * img_size) + (r * len) + c;
			t = (im * num_tomo * img_size) + (3 * img_size) + i;
			tomo[t] = img[p];
		}
/*				
		//Rot0, Flip1
		for(i = 0, r = 0; r < len; ++r)
		for(c = len - 1; c >= 0; --c, ++i)
		{
			p = (im * img_size) + (r * len) + c;
			t = (im * num_tomo * img_size) + (4 * img_size) + i;
			tomo[t] = img[p];
		}
		 		
		//Rot90, Flip1
		for(i = 0, c = 0; c < len; ++c)
		for(r = 0; r < len; ++r, ++i)
		{
			p = (im * img_size) + (r * len) + c;
			t = (im * num_tomo * img_size) + (5 * img_size) + i;
			tomo[t] = img[p];
		}

		//Rot180, Flip1
		for(i = 0, r = len - 1; r >= 0; --r)
		for(c = 0; c < len; ++c, ++i)
		{
			p = (im * img_size) + (r * len) + c;
			t = (im * num_tomo * img_size) + (6 * img_size) + i;
			tomo[t] = img[p];
		} 

		//Rot270, Flip1
		for(i = 0, c = len - 1; c >= 0; --c)
		for(r = len - 1; r >= 0; --r, ++i)
		{
			p = (im * img_size) + (r * len) + c;
			t = (im * num_tomo * img_size) + (7 * img_size) + i;
			tomo[t] = img[p];
		}
*/
	}
}
	
	
void free_mem() 
{
	free(img);
	free(single_photons), free(multi_photons);		
	free(tomo);
	free(noise);		
	free(avg_noise);
}
	
int poisson_dist( double m )
{
	int i = 0;
	double p, q, r;
		
	r = exp(-m);
	p = r;
	q = ((double) rand()) / RAND_MAX;
	
	while (p <= q)
	{
		++i;
		r *= m/i;
		p += r;
	}
	
	return i ;
}
	
double gaussian_dist(double mu, double var) 
{
	double p, q;
	do
	{
		p = ((double) rand()) / RAND_MAX;
		q = ((double) rand()) / RAND_MAX;	
	} while( (p < 1.E-14) && (q < 1.E-14) );
	
	return sqrt( - 2.0 * var * log( p ) ) * cos( 2 * PI * q ) + mu ;	
}


double exp_dist(double mu) 
{
	double p;
	p = ((double) rand()) / RAND_MAX;
	return (-1.*mu) * log(1 - p);
}
	
double fac(int n)
{
	int i;
	double s = 1.;
	if(n == 0)
		return s;
	else if (n > 0)
	{
		for(i = 1 ; i <= n ; i++)
			s *= i;
		return s;
	}
	else
	{
		fprintf(stderr, "Cannot compute negative factorials!\n");
		exit(1);
	}
}
