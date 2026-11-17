#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <stdbool.h>
#define Abs(a) ((a)>=0?(a):-(a))
#define R0 3   // mm
#define L0 25  // mm
#define Rf 2.45   // mm
#define Lf 28.78  // mm
#define GOODNESS 0.999899

static int total_observation(FILE *p) {
    int n=0;
    char line[50];
    while (fgets(line,50,p)!=NULL) {
        n++;
    }
    rewind(p);
    return n;
}

struct iteration {
    int *old_index;
    double *s;
    double *r;
};

static int compare_value(void const *a, void const *b) {
    struct iteration *it1= (struct iteration*) a;
    struct iteration *it2= (struct iteration*) b;

    if (*(it1->r) < *(it2->r))
        return 1;
    if (*(it1->r) > *(it2->r))
        return -1;
    return 0;
}

static int total_point(FILE *fp) {
    int n=0;
    char line[100];
    while (fgets(line,sizeof(line),fp)!=NULL) {
        n++;
    }
    rewind(fp);
    return n;
}

static double sum_x_pow(const double *x,int n,int p) {
    double s=0;
    for (int i=0; i<n; i++) {
        s+=pow(x[i],p*1.0);
    }
    return s;
}

static double sum_y_x_pow(const double *x,const double *y,const int n,const int p) {
    double s=0;
    for (int i=0; i<n; i++) {
        s+=pow(x[i],p*1.0)*y[i];
    }
    return s;
}

static double f(const double *co_eff, const double x) {
    return (co_eff[0]+ co_eff[1]*x+ co_eff[2]*x*x+ co_eff[3]*x*x*x+ co_eff[4]*x*x*x*x);
}

static double integrate(const double *co_eff, const double l, const double u) {
    int n= (u-l)/1e-6;
    double s=0;
    if (n%2==0) {
        s+=f(co_eff,l)+f(co_eff,u);
        for (int i=1; i<n; i++) {
            if ((i+1)%2==0) {
                s+=4*f(co_eff,l+i*1e-6);
            }
            else {
                s+=2*f(co_eff,l+i*1e-6);
            }
        }
        s=1e-6/3*s;
    }
    else {
        s+=3*1e-6/8*(f(co_eff,l)+ 3*f(co_eff,l+1e-6)+
            3*f(co_eff,l+2*1e-6)+ f(co_eff,l+3*1e-6));
        double s2=f(co_eff,l+3*1e-6)+f(co_eff,u);
        for (int i=4; i<n; i++) {
            if ((i+1)%2==0) {
                s2+=4*f(co_eff,l+i*1e-6);
            }
            else {
                s2+=2*f(co_eff,l+i*1e-6);
            }
        }
        s=s+1e-6/3*s2;
    }
    return s;
}

static void solve_equation_without_file(const int n,double **matrix,double *soln) {
    double *scaled_vec= malloc(sizeof(double)*n);
    if (!scaled_vec) {
        return;
    }
    for (int i=0; i<n; i++) {
        double mxm=0;
        for (int j=0; j<n; j++) {
            if (mxm<Abs(matrix[i][j]))
                mxm=Abs(matrix[i][j]);
        }
        scaled_vec[i]=mxm;
    }

    for (int i=0; i<n; i++) {
        int it=0;
        double *ratio=malloc(sizeof(double)*(n-i));
        if (!ratio) {
            free(scaled_vec);
            return;
        }
        double *new_scaled=malloc(sizeof(double)*(n-i));
        if (!new_scaled) {
            free(scaled_vec);
            free(ratio);
            return;
        }

        for (int j=i; j<n; j++) {
            new_scaled[it++]=scaled_vec[j];
        }
        it=0;
        for (int j=i; j<n; j++) {
            ratio[it]=Abs(matrix[j][i]/new_scaled[it]);
            it++;
        }

        int *old=malloc(sizeof(int)*(n-i));
        if (!old) {
            free(scaled_vec);
            free(ratio);
            free(new_scaled);
            return;
        }
        for (int k=i; k<n; k++) {
            old[k-i]=k;
        }

        struct iteration *iter=malloc(sizeof(struct iteration)*(n-i));
        if (!iter) {
            free(scaled_vec);
            free(ratio);
            free(new_scaled);
            free(old);
            return;
        }

        for (int j=0; j<n-i; j++) {
            iter[j].r=&ratio[j];
            iter[j].s=&new_scaled[j];
            iter[j].old_index=&old[j];
        }

        qsort(iter,n-i,sizeof(struct iteration),compare_value);

        int new=*iter[0].old_index;
        double *temp=malloc(sizeof(double)*(n+1));
        if (!temp) {
            free(scaled_vec);
            free(ratio);
            free(new_scaled);
            free(old);
            free(iter);
            return;
        }
        for (int x=0; x<n+1; x++) {
            temp[x]=matrix[i][x];
        }
        for (int x=0; x<n+1; x++) {
            matrix[i][x]=matrix[new][x];
        }
        for (int x=0; x<n+1; x++) {
            matrix[new][x]=temp[x];
        }
        double val=matrix[i][i];
        if (Abs(val)<1e-9) {
            printf("Singular Matrix\n");
            free(scaled_vec);
            free(ratio);
            free(new_scaled);
            free(old);
            free(iter);
            free(temp);
            return;
        }
        for (int r=0; r<n+1; r++) {
            matrix[i][r]/=val;
        }

        double tmp=scaled_vec[new];
        scaled_vec[new]=scaled_vec[i];
        scaled_vec[i]=tmp;

        for (int k=0; k<n; k++) {
            if (i==k)
                continue;
            double factor= -matrix[k][i]/matrix[i][i];
            for (int col=0; col<n+1; col++) {
                matrix[k][col]+=factor*matrix[i][col];
            }
        }

        free(ratio);
        free(iter);
        free(old);
        free(temp);
        free(new_scaled);
    }
    for (int i=0; i<n; i++) {
        soln[i]=matrix[i][n];
    }
    free(scaled_vec);
}

static void polyfit_without_file(double *x, double *y,double *soln,
    const int n, const int fit_rank) {
    double **matrix= malloc(sizeof(double*)*(fit_rank+1));
    if (!matrix) {
        return;
    }
    for (int i=0; i<=fit_rank; i++) {
        matrix[i]=malloc(sizeof(double)*(fit_rank+2));
        if (!matrix[i]) {
            free(matrix);
            return;
        }
    }

    for (int i=0; i<=fit_rank; i++) {
        for (int j=0; j<=fit_rank; j++) {
            matrix[i][j]=sum_x_pow(x,n,i+j);
        }
    }

    for (int i=0; i<=fit_rank; i++) {
        matrix[i][fit_rank+1]=sum_y_x_pow(x,y,n,i);
    }

    solve_equation_without_file(fit_rank+1,matrix,soln);

    // freeing here
    for (int i=0; i<=fit_rank; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

static double curve_goodness(double *stress, double *strain, const int n) {
    double *co_eff=malloc(sizeof(double)*2);
    polyfit_without_file(strain,stress,co_eff,n,1);

    double s_stress=0;
    for (int i=0; i<n; i++) {
        s_stress+=stress[i];
    }
    const double avg_stress= s_stress/n;

    double *fitted_stress=malloc(sizeof(double)*n);
    for (int i=0; i<n; i++) {
        fitted_stress[i]= co_eff[1]*strain[i]+co_eff[0];
    }

    double st=0;
    for (int i=0; i<n; i++) {
        st+=(stress[i]-avg_stress)*(stress[i]-avg_stress);
    }

    double sr=0;
    for (int i=0; i<n; i++) {
        sr+=(stress[i]-fitted_stress[i])*(stress[i]-fitted_stress[i]);
    }

    return (1-sr/st);
}

static void fit_engineering_curve(FILE *fp,FILE *fo, double *e_co, double *p_co) {
    const int total_obs=total_observation(fp);
    double *stress,*strain;
    stress=malloc(sizeof(double)*total_obs);
    if (!stress) {
        return;
    }
    strain=malloc(sizeof(double)*total_obs);
    if (!strain) {
        free(stress);
        return;
    }

    for (int i=0; i<total_obs; i++) {
        double a,b;
        fscanf(fp,"%lf\t%lf",&a,&b);
        stress[i]=a/(M_PI*R0*R0);
        strain[i]=b/L0;
    }

    for (int i=0; i<total_obs; i++) {
        if (i!=total_obs-1)
            fprintf(fo,"%.6lf\t%.6lf\n",stress[i],strain[i]);
        else
            fprintf(fo,"%.6lf\t%.6lf",stress[i],strain[i]);
    }

    int proportional_limit=0;

    for (int i=1; ; i++) {
        if (curve_goodness(stress,strain,total_obs*i/100)<=GOODNESS) {
            proportional_limit=total_obs*(i-0.5)/100;
            break;
        }
    }

    double proportional_stress=stress[proportional_limit-1];
    double proportional_strain=strain[proportional_limit-1];

    polyfit_without_file(strain,stress,e_co,proportional_limit,1);

    double E=e_co[1];
    // printf("%.6lf\n",e_co[0]);

    double *plastic_stress, *plastic_strain;
    plastic_stress=malloc(sizeof(double)*(total_obs-proportional_limit+1));
    if (!plastic_stress) {
        free(stress);
        free(strain);
        return;
    }
    plastic_strain=malloc(sizeof(double)*(total_obs-proportional_limit+1));
    if (!plastic_strain) {
        free(stress);
        free(strain);
        free(plastic_stress);
        return;
    }

    for (int i=proportional_limit-1; i<total_obs; i++) {
        plastic_stress[i-proportional_limit+1]=stress[i];
        plastic_strain[i-proportional_limit+1]=strain[i];
    }

    polyfit_without_file(plastic_strain,plastic_stress,p_co,
        total_obs-proportional_limit+1,4);

    double l=proportional_strain;
    double u=strain[total_obs-1];
    double diff=(u-l)/(1e7-1);

    double yield_strength=0;
    for (int i=0; i<1e7; i++) {
        if ((E*(l+i*diff-0.002))>
            (p_co[0]+
                p_co[1]*(l+i*diff)+
                    p_co[2]*(l+i*diff)*(l+i*diff)+
                        p_co[3]*(l+i*diff)*(l+i*diff)*(l+i*diff)+
                           p_co[4]*(l+i*diff)*(l+i*diff)*(l+i*diff)*(l+i*diff))) {
            yield_strength=E*(l+i*diff-0.002);
            break;
        }
    }
    // printf("%.6lf\n",yield_strength);

    double uts=-10;
    for (int i=0; i<1e7; i++) {
        double val=p_co[0]+
                p_co[1]*(l+i*diff)+
                    p_co[2]*(l+i*diff)*(l+i*diff)+
                        p_co[3]*(l+i*diff)*(l+i*diff)*(l+i*diff)+
                           p_co[4]*(l+i*diff)*(l+i*diff)*(l+i*diff)*(l+i*diff);
        if (uts<val) {
            uts=val;
        }
    }
    // printf("%.6lf\n",uts);

    printf("Proportional Limit stress : %.6lf\n",proportional_stress);
    printf("Proportional Limit strain : %.6lf\n",proportional_strain);
    printf("Young's Modulus           : %.6lf\n",E);
    printf("Yield Strength            : %.6lf\n",yield_strength);
    printf("Ultimate Tensile Strength : %.6lf\n",uts);
    printf("Fracture stress           : %.6lf\n",stress[total_obs-1]);
    printf("Fracture strain           : %.6lf\n",strain[total_obs-1]);
    printf("Modulus of Resilience     : %.6lf\n",0.5*proportional_strain*proportional_stress);
    printf("Modulus of Toughness      : %.6lf\n",0.5*proportional_strain*proportional_stress+
        integrate(p_co,l,u));

    free(stress);
    free(strain);
    free(plastic_stress);
    free(plastic_strain);
}

static double sum_xy(double *x, double *y, const int n) {
    double s=0;
    for (int i=0; i<n; i++) {
        s+=x[i]*y[i];
    }
    return s;
}

static double sum_x2(double *x, const int n) {
    double s=0;
    for (int i=0; i<n; i++) {
        s+=x[i]*x[i];
    }
    return s;
}

static void fit_true_curve(FILE *fp, FILE *fo, double *e_co, double *p_co,
    double *k, double *n, double *true_co_eff) {
    const int total_obs=total_observation(fp);
    double *stress,*strain;
    stress=malloc(sizeof(double)*total_obs);
    if (!stress) {
        return;
    }
    strain=malloc(sizeof(double)*total_obs);
    if (!strain) {
        free(stress);
        return;
    }

    for (int i=0; i<total_obs; i++) {
        double a,b;
        fscanf(fp,"%lf\t%lf",&a,&b);
        stress[i]=a/(M_PI*R0*R0);
        strain[i]=b/L0;
    }

    int proportional_limit=0;
    for (int i=1; ; i++) {
        if (curve_goodness(stress,strain,total_obs*i/100)<=GOODNESS) {
            proportional_limit=total_obs*(i-0.5)/100;
            break;
        }
    }

    // double proportional_stress=stress[proportional_limit-1];
    double proportional_strain=strain[proportional_limit-1];

    polyfit_without_file(strain,stress,e_co,proportional_limit,1);

    double E=e_co[1];
    // printf("%.6lf\n",e_co[0]);

    double *plastic_stress, *plastic_strain;
    plastic_stress=malloc(sizeof(double)*(total_obs-proportional_limit+1));
    if (!plastic_stress) {
        free(stress);
        free(strain);
        return;
    }
    plastic_strain=malloc(sizeof(double)*(total_obs-proportional_limit+1));
    if (!plastic_strain) {
        free(stress);
        free(strain);
        free(plastic_stress);
        return;
    }

    for (int i=proportional_limit-1; i<total_obs; i++) {
        plastic_stress[i-proportional_limit+1]=stress[i];
        plastic_strain[i-proportional_limit+1]=strain[i];
    }

    polyfit_without_file(plastic_strain,plastic_stress,p_co,
        total_obs-proportional_limit+1,4);

    double l=proportional_strain;
    double u=strain[total_obs-1];
    double diff=(u-l)/(1e7-1);

    double yield_strength=0;
    double yield_strain=0;
    for (int i=0; i<1e7; i++) {
        if ((E*(l+i*diff-0.002))>
            (p_co[0]+
                p_co[1]*(l+i*diff)+
                    p_co[2]*(l+i*diff)*(l+i*diff)+
                        p_co[3]*(l+i*diff)*(l+i*diff)*(l+i*diff)+
                           p_co[4]*(l+i*diff)*(l+i*diff)*(l+i*diff)*(l+i*diff))) {
            yield_strength=E*(l+i*diff-0.002);
            yield_strain=l+i*diff;
            break;
        }
    }
    // printf("%.6lf\n",yield_strength);

    double uts=-10;
    double uts_strain=0;
    for (int i=0; i<1e7; i++) {
        double val=p_co[0]+
                p_co[1]*(l+i*diff)+
                    p_co[2]*(l+i*diff)*(l+i*diff)+
                        p_co[3]*(l+i*diff)*(l+i*diff)*(l+i*diff)+
                           p_co[4]*(l+i*diff)*(l+i*diff)*(l+i*diff)*(l+i*diff);
        if (uts<val) {
            uts=val;
            uts_strain=l+i*diff;
        }
    }
    // printf("%.6lf\n",uts);

    int cnt=0;
    for (int i=0; i<total_obs; i++) {
        if (strain[i]<=uts_strain) {
            cnt++;
        }
    }

    double *true_stress= malloc(sizeof(double)*(cnt+1));
    if (!true_stress) {
        free(stress);
        free(strain);
        free(plastic_stress);
        free(plastic_strain);
        return;
    }

    double *true_strain= malloc(sizeof(double)*(cnt+1));
    if (!true_strain) {
        free(stress);
        free(strain);
        free(plastic_stress);
        free(plastic_strain);
        free(true_stress);
        return;
    }

    for (int i=0; i<cnt; i++) {
        true_stress[i]=stress[i]*(1+strain[i]);
        true_strain[i]=log(1+strain[i]);
    }

    true_stress[cnt]=stress[total_obs-1]*(R0*R0)/(Rf*Rf);
    true_strain[cnt]=(Lf-L0)/L0;

    for (int i=0; i<cnt; i++) {
        if (i!=cnt-1)
            fprintf(fo,"%.6lf\t%.6lf\n",true_stress[i],true_strain[i]);
        else
            fprintf(fo,"%.6lf\t%.6lf",true_stress[i],true_strain[i]);
    }

    double true_E= sum_xy(true_strain,true_stress,proportional_limit)/
        sum_x2(true_strain,proportional_limit);

    double true_yield= yield_strength*(1+yield_strain);
    double true_uts= uts*(1+uts_strain);
    double true_uts_strain= log(1+uts_strain);

    double *x= malloc(sizeof(double)*(cnt-proportional_limit));
    if (!x) {
        free(stress);
        free(strain);
        free(plastic_stress);
        free(plastic_strain);
        free(true_stress);
        free(true_strain);
        return;
    }
    double *y= malloc(sizeof(double)*(cnt-proportional_limit));
    if (!y) {
        free(stress);
        free(strain);
        free(plastic_stress);
        free(plastic_strain);
        free(true_stress);
        free(true_strain);
        free(x);
        return;
    }

    for (int i=proportional_limit; i<cnt; i++) {
        x[i-proportional_limit]=log(true_strain[i]);
        y[i-proportional_limit]=log(true_stress[i]);
    }

    double *eff= malloc(sizeof(double)*2);
    if (!eff) {
        free(stress);
        free(strain);
        free(plastic_stress);
        free(plastic_strain);
        free(true_stress);
        free(true_strain);
        free(x);
        return;
    }

    polyfit_without_file(x,y,eff,cnt-proportional_limit,1);
    free(x);
    free(y);

    *k=pow(M_E,eff[0]);
    *n=eff[1];

    if (proportional_limit>5) {
        x= malloc(sizeof(double)*5);
        y= malloc(sizeof(double)*5);
        for (int i=cnt-4; i<=cnt; i++) {
            x[i-cnt+4]= true_strain[i];
            y[i-cnt+4]= true_stress[i];
        }
        polyfit_without_file(x,y,true_co_eff,5,2);
    }
    else {
        x= malloc(sizeof(double)*2);
        y= malloc(sizeof(double)*2);
        for (int i=cnt-1; i<=cnt; i++) {
            x[i-cnt+1]= true_strain[i];
            y[i-cnt+1]= true_stress[i];
        }
        polyfit_without_file(x,y,true_co_eff,2,2);
    }

    printf("K = %.6lf\n",*k);
    printf("n = %.6lf\n",*n);
    printf("True Young's Modulus           : %.6lf\n",true_E);
    printf("True Ultimate Tensile Strength : %.6lf\n",true_uts);
    printf("True Ultimate Tensile Strain   : %.6lf\n",true_uts_strain);
    printf("True Fracture Stress           : %.6lf\n",true_stress[cnt]);
    printf("True Fracture Strain           : %.6lf\n",true_strain[cnt]);

    free(stress);
    free(strain);
    free(plastic_stress);
    free(plastic_strain);
    free(true_stress);
    free(true_strain);
    free(eff);
    free(x);
    free(y);
}

int main() {
    FILE *fp=fopen("double_yield_point_tensile_data.txt","r");
    if (!fp) {
        return 1;
    }

    // FILE *fo=fopen("single_yield_point_stress_strain.txt","w");
    // if (!fo) {
    //     fclose(fp);
    //     return 1;
    // }
    //
    double *elastic_eff= malloc(sizeof(double)*2);
    if (!elastic_eff) {
        fclose(fp);
        // fclose(fo);
        return 1;
    }

    double *plastic_eff= malloc(sizeof(double)*5);
    if (!plastic_eff) {
        fclose(fp);
        // fclose(fo);
        free(elastic_eff);
        return 1;
    }

    // fit_engineering_curve(fp,fo,elastic_eff,plastic_eff);
    // printf("%.6lf\n",elastic_eff[1]);
    // for (int i=0; i<5; i++) {
    //     printf("%.6lf\n",plastic_eff[i]);
    // }

    FILE *fo= fopen("double_yield_point_true_data.txt","w");
    if (!fo) {
        fclose(fp);
        return 1;
    }

    double *true_co_eff= malloc(sizeof(double)*3);
    if (!true_co_eff) {
        fclose(fp);
        fclose(fo);
        return 1;
    }

    double *k= malloc(sizeof(double));
    double *n= malloc(sizeof(double));
    fit_true_curve(fp,fo,elastic_eff,plastic_eff,k,n,true_co_eff);
    for (int i=0; i<3; i++) {
        printf("%.6lf\n",true_co_eff[i]);
    }

    fclose(fp);
    fclose(fo);
    free(true_co_eff);
    free(k);
    free(n);
    free(elastic_eff);
    free(plastic_eff);
}