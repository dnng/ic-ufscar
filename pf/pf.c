#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/*
 * Valid numbers: 100, 500, 2000, and 10000.
 * Maximum is 10000. Doesn't convege with 1000.
 */
#define MAX_PARTICLES 2000

/*
 * Reads sensor data from file. MUST BE DEFINED!!!
 */
#define BENCH_DATA

#define PI 3.14159265

#if 0
/*
 * Replaces rand() function by values stored in a file.
 */
#define BENCH_RAND
#endif

/*
 * Maps main structure
 */
typedef struct st_map MAP;
struct st_map
{
    int size_x;
    int size_y;
    int grid[561][208];
    int lh_grid[561][208];
    float scale;
};

/*
 * Particle structure
 */
typedef struct st_part PARTICLE;
struct st_part
{
    float x, y, a;
    unsigned w;
    unsigned accw;
    int valid;
};

MAP map;
PARTICLE particle[MAX_PARTICLES];

/*
 * replaces rand function when BENCH_RAND defined
 */
int bench_rand(FILE *rand_f)
{
    int r;
#ifdef BENCH_RAND
    fscanf(rand_f, "%i ", &r);
#else
    r = rand();
#endif
    return r;
}

/*
 * load laser angles from file
 */
void bench_read_laser_ang( FILE *laser_ang_f, double laser[181][2])
{
    int i;
    float ang;

    for(i=0; i<181; i++) {
        fscanf(laser_ang_f, "%f", &ang);
        laser[i][1] = ang;
    }
}

/*
 * reads sensor data from file. Must define BENCH_DATA
 */
void bench_read_data(FILE *position_f, FILE *laser_f,
                     double position[3], double laser[181][2])

{
    int i;
    float px, py, pa, scan;

    fscanf(position_f, "%f %f %f\n", &px, &py, &pa);

    position[0] = px;
    position[1] = py;
    position[2] = pa;

    for(i = 0; i < 181; i++)
    {
        fscanf(laser_f, "%f", &scan);
        laser[i][0] = scan;
    }
}

/*
 * round angles
 */
inline float round_angle(float a)
{
    while (a >  PI) {
        a -= 6.2831853;
    }
    while (a < -PI) {
        a += 6.2831853;
    }
    return a;
}

/*
 * sample from a normal function, given standard deviation
 */
float normal(float sd, FILE *rand_f)
{
    int i;
    float num;

    num = 0;
    for(i = 0;i < 12;i++) {
        num += (bench_rand(rand_f) / (RAND_MAX/2.0)) - 1;
    }

    num = num * sd / 6;
    return num;
}

/*
 * writes particles to a file (particles_it[number].txt)
 */
void store_particles(int iteration)
{
    register int i;
    FILE * stream;
    char particle_file[20];

    sprintf(particle_file, "%s%d.txt", "particle_it",iteration);
    stream = fopen(particle_file, "w");

#if 0
    int coord_1, coord_2;
    int invalid = 0;
    for (i=0;i<MAX_PARTICLES;i++) {
        coord_1 = (int) (particle[i].x*map.scale);
        coord_2 = (int) (map.size_y-(particle[i].y*map.scale))
        int res = map.grid[coord_1][coord_2];
        if (res != -1) {
            particle[i].valid=0;
            invalid++;
        } else {
            particle[i].valid = 1;
        }
        if (particle[i].valid){
            fprintf(stream, "%f %f\n", particle[i].x, particle[i].y);
        }
    }
#endif

    for (i = 0;i < MAX_PARTICLES;i++){
        fprintf(stream, "%f %f\n", particle[i].x, particle[i].y);
    }
}

#if 0
/*
 * calc x, y, and a mean and stdev
 *
 * Commented because of the test of the store_particles() function
 */
void calc_mean() {
    register int i;
    double mean_x  = 0;
    double mean_y  = 0;
    double stdev_x = 0;
    double stdev_y = 0;
    int invalid    = 0;
    int coord_1    = 0;
    int coord_2    = 0;
    int res        = 0;

    double mean_sin = 0;
    double mean_cos = 0;
    double mean_a;

    for (i = 0;i < MAX_PARTICLES;i++)
    {
        coord_1 = (int) (particle[i].x*map.scale);
        coord_2 = (int) (map.size_y-(particle[i].y*map.scale));
        res = map.grid[coord_1][coord_2];
        if (res != -1) {
            particle[i].valid= 0;
            invalid++;
        }
        else {
            particle[i].valid = 1;
        }
    }

    for (i=0;i<MAX_PARTICLES;i++) {
        if (particle[i].valid) {
            mean_x    += particle[i].x;
            mean_y    += particle[i].y;
            mean_sin  += sin(particle[i].a);
            mean_cos  += cos(particle[i].a);
        }
    }

    mean_x   /= (MAX_PARTICLES - invalid);
    mean_y   /= (MAX_PARTICLES - invalid);
    mean_sin /= (MAX_PARTICLES - invalid);
    mean_cos /= (MAX_PARTICLES - invalid);

    mean_a = atan2(mean_sin, mean_cos);

    for (i = 0; i < MAX_PARTICLES; i++) {
        if (particle[i].valid) {
            stdev_x += fabsf(particle[i].x - mean_x);
            stdev_y += fabsf(particle[i].y - mean_y);
        }
    }

    stdev_x /= (MAX_PARTICLES - invalid);
    stdev_y /= (MAX_PARTICLES - invalid);

    printf("\rMean x:%.2lf , y:%.2lf , a:%.2lf StdDev x:%.5lf, y:%.5lf",
            mean_x,  mean_y, mean_a*57.3,  stdev_x, stdev_y);
}
#endif

/* load map from file
 */
void load_map(const char *file_name)
{
    int i, j, aux;
    FILE * file;

    file = fopen("rth.map", "r");
    if (file != NULL) {
        fscanf(file, "mapd %i %i", &map.size_x, &map.size_y);
        printf("\nLoading map with dimensions x:%i, y:%i", map.size_x, map.size_y);

        for(i = 0; i < map.size_y; i++) {
            /*
             * x and y inverted
             */
            for(j = 0; j < map.size_x; j++) {
                fscanf(file,"%i\n", &aux);
                switch(aux) {
                    case 1:
                        map.grid[j][i] = 1;
                        break;
                    case 0:
                        map.grid[j][i] = -1;
                        break;
                    case -1:
                        map.grid[j][i] = 0;
                        break;
                }
            }
        }
    }
    else {
        printf("\nFile Error\nTrying to open:'%s'\n",file_name);
    }
}


/*
 * calculates likelyhood map
 */
void calc_lh_map(void)
{
    int i,j,k,dist;
    int aux1s, aux1f, aux2s, aux2f;
    int lh_value[] = {50, 7, 1, 0, 0, 0, 0};

    printf("\nClaculating Likelyhood Map...");

    for(i = 0; i < map.size_x; i++) {
        for(j = 0; j < map.size_y; j++) {
            map.lh_grid[i][j] = 0;
        }
    }


    for(dist = 6; dist >= 0; dist--) {
        for(i = 0; i < map.size_x; i++) {
            for(j = 0; j < map.size_y; j++) {
                aux1s = i - dist;
                aux1f = i + dist;
                aux2s = j - dist;
                aux2f = j + dist;

                if (aux1s < 0)
                    aux1s = 0;

                if (aux2s < 0)
                    aux2s = 0;

                if (aux1f >= map.size_x)
                    aux1f = map.size_x - 1;

                if (aux2f >= map.size_y)
                    aux2f = map.size_y - 1;

                for (k = aux1s; k <= aux1f; k++) {
                    if (map.grid[k][aux2s] == 1)
                        map.lh_grid[i][j] = lh_value[dist];
                    if (map.grid[k][aux2f] == 1)
                        map.lh_grid[i][j] = lh_value[dist];
                }

                for (k = aux2s; k <= aux2f; k++) {
                    if (map.grid[aux1s][k] == 1)
                        map.lh_grid[i][j] = lh_value[dist];
                    if (map.grid[aux1f][k] == 1)
                        map.lh_grid[i][j] = lh_value[dist];
                }
            }
        }
    }
    printf(" Done!\n");
}

/*returns a random particle position */
inline void rand_position(float *rx, float *ry, float *ra, FILE *rand_f)
{
    int coord_1 = 0;
    int coord_2 = 0;
    *rx = (bench_rand(rand_f) % map.size_x*100) / (map.scale*100);
    *ry = (bench_rand(rand_f) % map.size_y*100) / (map.scale*100);
    *ra = ((bench_rand(rand_f)%62)-31) / 10.0;

    coord_1 = (int)((*rx) * map.scale);
    coord_2 = map.size_y- ((int)((*ry) * map.scale));

    while (map.grid[coord_1][coord_2] != -1)
    {
        *rx = (bench_rand(rand_f) % map.size_x*100)/(map.scale*100);
        *ry = (bench_rand(rand_f) % map.size_y*100)/(map.scale*100);

        coord_1 = (int)((*rx) * map.scale);
        coord_2 = map.size_y- ((int)((*ry) * map.scale));
    }
}

/*
 * resample one particle
 */
inline void resample_one(int i, FILE *rand_f)
{
    int ref;
    int coord_1 = 0;
    int coord_2 = 0;

    ref     = bench_rand(rand_f) % MAX_PARTICLES;
    coord_1 = (int) (particle[ref].x * map.scale);
    coord_2 = (int) (map.size_y - (particle[ref].y * map.scale));

    while (map.grid[coord_1][coord_2] != -1)
    {
        ref = bench_rand(rand_f) % MAX_PARTICLES;
        coord_1 = (int) (particle[ref].x * map.scale);
        coord_2 = (int) (map.size_y - (particle[ref].y * map.scale));
    }

    particle[i].x = particle[ref].x;
    particle[i].y = particle[ref].y;
    particle[i].a = particle[ref].a;
}

/*
 * resample all particles
 */
inline void resample_all(unsigned total, float res_rate, FILE *rand_f)
{
    int i;
    long r, j, max_j, min_j, old_j;
    float rx, ry, ra;

    for(i = 0; i < (int) MAX_PARTICLES * res_rate; i++) {
        r     = (long) (bench_rand(rand_f)%(total));
        j     = MAX_PARTICLES/2;
        max_j = (MAX_PARTICLES-1);
        min_j = 0;

        while (1) {
            if ((j > 0) && (j < (MAX_PARTICLES - 1))) {
                if ((particle[j-1].accw < r) && (particle[j].accw >= r))
                    break;
            }
            else
                break;

            if (particle[j].accw < r) {
                min_j = j;
                old_j = j;
                j     = (min_j+max_j)/2;
                if (j == old_j) {
                    min_j++;
                    j++;
                }
            }
            else {
                max_j = j;
                old_j = j;
                j     = (min_j + max_j) / 2;
                if (j == old_j) {
                    max_j--;
                    j--;
                }
            }
        }

        particle[i].x = particle[j].x;
        particle[i].y = particle[j].y;
        particle[i].a = particle[j].a;
    }

    for(i = (int) MAX_PARTICLES*res_rate; i < MAX_PARTICLES; i++) {
        rand_position(&rx, &ry, &ra, rand_f);
        particle[i].x = rx;
        particle[i].y = ry;
        particle[i].a = ra;
    }
}

/*
 * Init particles
 */
void init_particles2(float x, float y, float a) {
    printf("\nInitializing particles...");
    int i;

#ifdef BENCH_DATA
    FILE *part_init;
    float init_x, init_y, init_a;

    part_init = fopen("part_init.dat", "r");

    for (i = 0; i < MAX_PARTICLES; i++) {
        fscanf(part_init,"%f %f %f\n", &init_x, &init_y, &init_a);
        particle[i].x = init_x;
        particle[i].y = init_y;
        particle[i].a = init_a;
        particle[i].w = 0;
    }
    fclose(part_init);
#else
    float var_xy=8.0, var_a=0.8;

    for (i=0; i<MAX_PARTICLES; i++) {
        particle[i].x = x + (((rand()%2000*var_xy)/1000.0)-var_xy);
        particle[i].y = y + (((rand()%2000*var_xy)/1000.0)-var_xy);
        particle[i].a = a + (((rand()%2000*var_a) /1000.0)- var_a);
        particle[i].w = 0;
    }
#endif
    printf("Done.");
}

/*
 * Calculates the observation model (weight particles)
 */
inline unsigned do_observation_model(float x, float y, float a, double laser[181][2])
{
    int i, res = 0;
    int gx, gy;

    gx = (int) x*map.scale;
    gy = (int) map.size_y-(y*map.scale);

    if ((gx<0) || (gx>=map.size_x) || (gy<0) || (gy>=map.size_y)) {
        return 0;
    }

    if (map.grid[gx][gy] != -1) {
        return 0;
    }

    for(i = 0; i < 181; i += 10) {
        gx = (x + laser[i][0] * cos(a + laser[i][1])) * map.scale;
        gy = map.size_y-((y+laser[i][0]*sin(a+laser[i][1])) * map.scale);

        if ((gx<0) || (gx>=map.size_x) || (gy<0) || (gy>=map.size_y))
            continue;

        if ((laser[i][0]>7.98) && (map.grid[gx][gy]==-1))
            res+=2;
        else
            res+=map.lh_grid[gx][gy];
    }

    return res/10;
}

/*
 * Perform the action model (propagate particles)
 */
inline void do_action_model3(float odom_ang, float odom_tr, float xo,
                             float yo, float to, float *xr,
                             float *yr, float *tr, FILE *rand_f)
{
    float odom_ang_b, odom_tr_b;
    float gain_t = 0.8, gain_r = 0.5;
    float norm_t, norm_r;

    norm_t = normal(0.8, rand_f);
    norm_r = normal(0.7, rand_f);

    odom_tr_b =
        odom_tr +
        (odom_tr * norm_t * gain_t) +
        (odom_ang * norm_r * gain_r * 0.5);

    odom_ang_b =
        odom_ang +
        (odom_ang * norm_r * gain_r) +
        (odom_tr  * norm_r * gain_r);

    *tr = to + odom_ang_b;
    *tr = round_angle(*tr);

    *xr = xo + (odom_tr_b*cos(*tr));
    *yr = yo + (odom_tr_b*sin(*tr));
}


/* Main loop */
void filter_loop(void)
{
    float xo, yo, ao;
    float dif_x, dif_y, dif_a;
    int i, count=0, invalid;
    unsigned obs_total;
    float odom_ang, odom_tr;

    double position[3];
    double laser[181][2];

    int skip_obs  =5;

    FILE *position_f, *laser_f, *laser_ang_f;
    FILE *rand_f;

    position_f=fopen("position.dat", "r");
    laser_f=fopen("laser.dat", "r");
    laser_ang_f=fopen("laser_ang.dat", "r");

#ifdef BENCH_RAND
    rand_f = fopen("rand.dat", "r");
#else
    rand_f = NULL;
#endif
    printf("\nEntering loop\n");
    printf("\nCreating files:\n");
#ifdef BENCH_DATA
    bench_read_laser_ang(laser_ang_f, laser);
    bench_read_data(position_f, laser_f, position, laser);
#else

    for(i=0; i<5; i++)
        playerc_client_read(client);
#endif
    xo = position[0];
    yo = position[1];
    ao = position[2];

    int it=0;
    /*here was count<200 instead of count<1 */
    while(count<200)
    {
        /*
         * Commented because of the test
         * of the store_particles function
         *
         calc_mean();
         */
        if((count%10) == 0){
            it++;
            store_particles(it);
            printf("File %d\n", it);
        }
        count++;

#ifdef BENCH_DATA
        bench_read_data(position_f, laser_f, position, laser);
#else
        playerc_client_read(client);
        playerc_client_read(client);
        playerc_client_read(client);
#endif
        dif_x = position[0] - xo;
        dif_y = position[1] - yo;
        dif_a = position[2] - ao;
        dif_a=round_angle(dif_a);

        if ((int)(dif_x*1000) || (int)(dif_y*1000) || (int)(dif_a*1000)) {
            odom_tr  = sqrt((dif_x*dif_x)+(dif_y*dif_y));
            odom_ang = dif_a;
            odom_ang = round_angle(odom_ang);

            xo = position[0];
            yo = position[1];
            ao = position[2];

            invalid=0;

            for (i=0; i < MAX_PARTICLES; i++) {

                do_action_model3(odom_ang, odom_tr, particle[i].x,
                        particle[i].y, particle[i].a, &particle[i].x,
                        &particle[i].y, &particle[i].a, rand_f);

                if (map.grid[(int) (particle[i].x*map.scale)]
                    [(int) (map.size_y-(particle[i].y*map.scale))] != -1) {
                    resample_one(i, rand_f);
                    invalid++;
                }
                particle[i].a = round_angle(particle[i].a);
            }

            if ((count % skip_obs)==0) {
                obs_total=0;

                for (i=0; i < MAX_PARTICLES; i++) {
                    particle[i].w = do_observation_model(particle[i].x,
                            particle[i].y, particle[i].a, laser);
                    obs_total += particle[i].w;
                    particle[i].accw = obs_total;
                }

                resample_all(obs_total, 1.0, rand_f);
            }
        }
    }

    fclose(position_f);
    fclose(laser_f);
    fclose(laser_ang_f);

#ifdef BENCH_RAND
    fclose(rand_f);
#endif
}

int main(void)
{
    const char * map_file = "rth.map";
    load_map(map_file);
    map.scale = 10;
    calc_lh_map();
    init_particles2(25.30, 3.10, 0);
    filter_loop();
    return 0;
}
