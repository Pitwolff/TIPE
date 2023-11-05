#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>
#include <time.h>


#define nb_params 10
#define pi 3.14159265358979323846

struct individu {double* data; int score;};
typedef struct individu ind;

struct mutateur {double* data; int nb_enfants; ind** enfants; int score;};
typedef struct mutateur mut;

struct population {ind** individus; int capacite;};
typedef struct population pop;

struct population2 {mut** individus; int capacite;};
typedef struct population2 pop2;

int max(int a, int b) {
    return (a > b) ? a : b;
}

double g1(ind* joueur) {
    double coord[nb_params] = {1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000};
    double d = 0;
    for (int i = 0; i != nb_params; i += 1) {
        d += (coord[i] - joueur->data[i]) * (coord[i] - joueur->data[i]);
    }
    assert(!isnan(d));
    return sqrt(d);
}

double f1(ind* j1, ind* j2) {
    return g1(j1) - g1(j2);
}

double pol(double* coeffs, double x) {
    double s = 0;
    for (int i = nb_params - 1; i != - 1; i -= 1) {
        s *= x;
        s += coeffs[i];
    }
    return s;
}

double g2(ind* joueur) {
    double coeffs[nb_params] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double d = 0;
    for (int i = 0; i != nb_params; i += 1) {
        d += pow((pol(coeffs, i) - pol(joueur->data, i)), 2);
    }
    assert(!isnan(d));
    return sqrt(d);
}

double f2(ind* j1, ind* j2) {
    return g2(j1) - g2(j2);
}

void tournoi(pop* joueurs, int nb_individus, double (*f) (ind*, ind*)) { // la taille doit être une puissance de 2 impaire et pas d'égalité
    int t = joueurs->capacite;
    int* classement1 = malloc(t * sizeof(int));
    int* classement2 = malloc(t * sizeof(int));
    for (int i = 0; i != t; i += 1) {
        classement1[i] = i;
    }
    for (int i = 0; i != t; i += 1) {
        joueurs->individus[i]->score = 0;
    }
    int p = 2 * nb_individus + 1;
    int* coeffs_binomiaux = malloc(p * sizeof(int));
    coeffs_binomiaux[0] = 1;
    for (int i = 1; i != p; i += 1) {
        coeffs_binomiaux[i] = 0;
    }
    for (int k = 0; k != p; k += 1) {
        for (int i = k; i != 0; i -= 1) {
            coeffs_binomiaux[i] = coeffs_binomiaux[i - 1] + coeffs_binomiaux[i];
        }
        int deb = 0;
        for (int j = 0; j != k + 1; j += 1) {
            int fin = deb + t * coeffs_binomiaux[j];
            int nb = (fin - deb) / 2;
            if (k % 2 == 0) {
                    for (int i = 0; i != nb; i += 1) {
                        if (f(joueurs->individus[classement1[deb + i]], joueurs->individus[classement1[fin - i - 1]]) < 0) {
                            joueurs->individus[classement1[deb + i]]->score += 1;
                            joueurs->individus[classement1[fin - i - 1]]->score -= 1;
                            classement2[deb + i] = classement1[deb + i];
                            classement2[deb + nb + i] = classement1[fin - i - 1];
                        }
                        else {
                            joueurs->individus[classement1[deb + i]]->score -= 1;
                            joueurs->individus[classement1[fin - i - 1]]->score += 1;
                            classement2[deb + i] = classement1[fin - i - 1];
                            classement2[deb + nb + i] = classement1[deb + i];
                        }
                    }
                }
            else {
                    for (int i = 0; i != nb; i += 1) {
                        if (f(joueurs->individus[classement2[deb + i]], joueurs->individus[classement2[fin - i - 1]]) < 0) {
                            joueurs->individus[classement2[deb + i]]->score += 1;
                            joueurs->individus[classement2[fin - i - 1]]->score -= 1;
                            classement1[deb + i] = classement2[deb + i];
                            classement1[deb + nb + i] = classement2[fin - i - 1];
                        }
                        else {
                            joueurs->individus[classement2[deb + i]]->score -= 1;
                            joueurs->individus[classement2[fin - i - 1]]->score += 1;
                            classement1[deb + i] = classement2[fin - i - 1];
                            classement1[deb + nb + i] = classement2[deb + i];
                        }
                    }
                } 
            deb = fin;
        }
        t /= 2;
    }
    free(classement1);
    free(classement2);
    free(coeffs_binomiaux);
}

void tournoi2(pop2* mutateurs) {
    for (int i = 0; i != mutateurs->capacite; i += 1) {
        mut* mutateur = mutateurs->individus[i];
        mutateur->score = INT_MIN;
        for (int j = 0; j != mutateurs->capacite; j += 1) {
            mutateur->score = max(mutateur->score, mutateur->enfants[j]->score);
        }
    }
}

double f_mutation() {
    double x = (double)rand() / RAND_MAX;
    while (x == 0 || x == 1) {
        x = (double)rand() / RAND_MAX;
    }
    return (log(x) - log(1 - x)) * sqrt(3) / pi;
}

ind* creer_ind() {
    double* data = malloc(sizeof(double) * nb_params);
    for (int i = 0; i != nb_params; i += 1) {
        data[i] = f_mutation();
    }
    ind* joueur = malloc(sizeof(ind));
    joueur->data = data;
    joueur->score = 0;
    return joueur;
}

mut* creer_mut(int nb_enfants) {
    double* data = malloc(sizeof(double) * nb_params);
    ind** enfants = malloc(sizeof(ind*) * nb_enfants);
    for (int i = 0; i != nb_params; i += 1) {
        data[i] = f_mutation();
    }
    mut* mutateur = malloc(sizeof(mut));
    mutateur->data = data;
    mutateur->nb_enfants = nb_enfants;
    mutateur->enfants = enfants;
    mutateur->score = INT_MIN;
    return mutateur;
}

void muter_ind(ind* joueur, mut* mutateur) {
    for (int i = 0; i != nb_params; i += 1) {
        joueur->data[i] += f_mutation() * mutateur->data[i];
    }
}

void muter_mut(mut* mutateur) {
    for (int i = 0; i != nb_params; i += 1) {
        mutateur->data[i] *= exp(f_mutation());
    }
}

ind* copy_ind(ind* joueur) {
    double* data = malloc(sizeof(double) * nb_params);
    for (int i = 0; i != nb_params; i += 1) {
        data[i] = joueur->data[i];
    }
    ind* new_joueur = malloc(sizeof(ind));
    new_joueur->data = data;
    new_joueur->score = 0;
    return new_joueur;
}

mut* copy_mut(mut* mutateur) {
    double* data = malloc(sizeof(double) * nb_params);
    ind** enfants = malloc(sizeof(ind*) * mutateur->nb_enfants);
    for (int i = 0; i != nb_params; i += 1) {
        data[i] = mutateur->data[i];
    }
    mut* new_mutateur = malloc(sizeof(mut));
    new_mutateur->data = data;
    new_mutateur->nb_enfants = mutateur->nb_enfants;
    new_mutateur->enfants = enfants;
    new_mutateur->score = INT_MIN;
    return new_mutateur;
}

void lib_mem_ind(ind* joueur) {
    free(joueur->data);
    free(joueur);
}

void lib_mem_mut(mut* mutateur) {
    free(mutateur->data);
    free(mutateur->enfants);
    free(mutateur);
}

pop* creer_pop(int capacite) {
    ind** individus = malloc(sizeof(ind*) * capacite);
    for (int i = 0; i != capacite; i += 1) {
        individus[i] = creer_ind();
    }
    pop* joueurs = malloc(sizeof(pop));
    joueurs->individus = individus;
    joueurs->capacite = capacite;
    return joueurs;
}

pop2* creer_pop2(int capacite) {
    mut** individus = malloc(sizeof(mut*) * capacite);
    for (int i = 0; i != capacite; i += 1) {
        individus[i] = creer_mut(capacite);
    }
    pop2* mutateurs = malloc(sizeof(pop2));
    mutateurs->individus = individus;
    mutateurs->capacite = capacite;
    return mutateurs;
}

void lib_mem_pop(pop* joueurs) {
    for (int i = 0; i != joueurs->capacite; i += 1) {
        lib_mem_ind(joueurs->individus[i]);
    }
    free(joueurs->individus);
    free(joueurs);
}

void lib_mem_pop2(pop2* mutateurs) {
    for (int i = 0; i != mutateurs->capacite; i += 1) {
        lib_mem_mut(mutateurs->individus[i]);
    }
    free(mutateurs->individus);
    free(mutateurs);
}

void selection(pop* joueurs) {
    int j = 0;
    for (int i = 0; i != joueurs->capacite; i += 1) {
        if (joueurs->individus[i]->score > 0) {
            joueurs->individus[j] = joueurs->individus[i];
            j += 1;
        }
        else {
            lib_mem_ind(joueurs->individus[i]);
        }
    }
}

void selection2(pop2* mutateurs, int nb_individus) {
    mut** new_individus = malloc(sizeof(mut*) * mutateurs->capacite);
    int c = 0;
    for (int i = 2 * nb_individus + 1; i != -2 * (nb_individus + 1); i -= 1) {
        for (int j = 0; j != mutateurs->capacite; j += 1) {
            if (mutateurs->individus[j]->score == i) {
                new_individus[c] = mutateurs->individus[j];
                c += 1;
            }
        }
    }
    for (int i = mutateurs->capacite / 2; i != mutateurs->capacite; i += 1) {
        lib_mem_mut(new_individus[i]);
    }
    free(mutateurs->individus);
    mutateurs->individus = new_individus;
}

void reproduction(pop* joueurs, pop2* mutateurs) {
    int nb = mutateurs->capacite;
    int nb2 = joueurs->capacite / 2;
    for (int i = 0; i != nb; i += 1) {
        for (int j = 0; j != nb; j += 1) {
            joueurs->individus[nb * i + j + nb2] = copy_ind(joueurs->individus[nb * i + j]);
            muter_ind(joueurs->individus[nb * i + j + nb2], mutateurs->individus[i]);
            mutateurs->individus[i]->enfants[j] = joueurs->individus[nb * i + j + nb2];
        }
    }
}

void reproduction2(pop2* mutateurs) {
    int nb = mutateurs->capacite / 2;
    for (int i = nb; i != mutateurs->capacite; i += 1) {
        mutateurs->individus[i] = copy_mut(mutateurs->individus[i - nb]);
        muter_mut(mutateurs->individus[i]);
    }
}

void evolution(pop* joueurs, pop2* mutateurs, int nb_individus, double (*f) (ind*, ind*)) {
    tournoi(joueurs, nb_individus, f);
    tournoi2(mutateurs);
    selection2(mutateurs, nb_individus);
    selection(joueurs);
    reproduction2(mutateurs);
    reproduction(joueurs, mutateurs);
}

ind* best_player(pop* joueurs) {
    double max = 0;
    int pos = -1;
    for (int i = 0; i != joueurs->capacite; i += 1) {
        if (joueurs->individus[i]->score > max) {
            max = joueurs->individus[i]->score;
            pos = i;
        }
    }
    return joueurs->individus[pos];
}

void afficher_joueur(ind* joueur) {
    for (int i = 0; i != nb_params; i += 1) {
        printf("%d %.10f\n", i, joueur->data[i]);
    }
    printf("perf : %.10f\n", g1(joueur)); // changer g1 en g2 si on utilise f2
}

mut* best_mutateur(pop2* mutateurs) {
    return mutateurs->individus[0];
}

void afficher_mutateur(mut* mutateur) {
    for (int i = 0; i != nb_params; i += 1) {
        printf("%d %.10f\n", i, mutateur->data[i]);
    }
    printf("score : %d\n", mutateur->score);
}

void afficher_couple(ind* joueur, mut* mutateur) {
    for (int i = 0; i != nb_params; i += 1) {
        printf("%d %.10f  %.10f\n", i, joueur->data[i], mutateur->data[i]);
    }
    printf("perf : %.10f  score : %d\n", g1(joueur), mutateur->score); // changer g1 en g2 si on utilise f2
}

ind* main2(double (*f) (ind*, ind*), int temps, int nb_individus) {
    int t0 = time(NULL);
    pop* joueurs = creer_pop(1<<(2 * nb_individus + 1));
    pop2* mutateurs = creer_pop2(1 << (nb_individus));
    for (int i = 0; i != mutateurs->capacite; i += 1) {
        for (int j = 0; j != mutateurs->capacite; j += 1) {
            mutateurs->individus[i]->enfants[j] = joueurs->individus[mutateurs->capacite * i + j];
        }
    }
    int i = 0;
    while(time(NULL) - t0 < temps) {
        evolution(joueurs, mutateurs, nb_individus, f);
        if (time(NULL) - t0 > temps * i / 10) {
            printf("%d\n", i);
            afficher_couple(best_player(joueurs), best_mutateur(mutateurs));
            i += 1;
        }
    }
    ind* meilleur_joueur = copy_ind(best_player(joueurs));
    lib_mem_pop(joueurs);
    lib_mem_pop2(mutateurs);
    return meilleur_joueur;
}

int main() {
    srand(time(NULL));
    ind* joueur = main2(f1, 1, 3);
    //ind* joueur = main2(f2, 10, 6);
    afficher_joueur(joueur);
    lib_mem_ind(joueur);
    return 0;
}
