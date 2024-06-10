#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>
#include <time.h>

#define pi 3.14159265358979323846

struct mutateur;

struct individu {double* data; int score; struct mutateur* pere; int nb_params;};
typedef struct individu ind;

struct mutateur {double* data; int nb_enfants; ind** enfants; int score; int nb_enfants_en_vie; int nb_params;};
typedef struct mutateur mut;

struct population2 {mut** individus; int capacite; int log2_capacite;};
typedef struct population2 pop2;

struct population {ind** individus; int capacite; int log2_capacite; pop2* mutateurs; int nb_vic;};
typedef struct population pop;

struct populations {pop** individus; int plus_faible; int* victoires; int nb_params;};
typedef struct populations pops;

int max(int a, int b) {
    return (a > b) ? a : b;
}

double rand_normal(double ecart_type) {
    double u1 = (rand() + 1.0) / (RAND_MAX + 2.0);
    double u2 = rand() / (RAND_MAX + 1.0);
    double z = sqrt(-2.0 * log(u1)) * cos(2.0 * pi * u2);
    return ecart_type * z;
}

void tournoi(pop* joueurs, bool (*f) (double*, double*)) { // la taille doit être une puissance de 2 et pas d'égalité
    int t = joueurs->capacite;
    int* classement1 = malloc(t * sizeof(int));
    int* classement2 = malloc(t * sizeof(int));
    for (int i = 0; i != t; i += 1) {
        classement1[i] = i;
    }
    for (int i = 0; i != t; i += 1) {
        joueurs->individus[i]->score = 0;
    }
    int p = joueurs->log2_capacite;
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
                        if (f(joueurs->individus[classement1[deb + i]]->data, joueurs->individus[classement1[fin - i - 1]]->data)) {
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
                        if (f(joueurs->individus[classement2[deb + i]]->data, joueurs->individus[classement2[fin - i - 1]]->data)) {
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
    ind** new_individus = malloc(sizeof(ind*) * joueurs->capacite);
    if (p && 1 == 0) {
        for (int i = 0; i != joueurs->capacite; i += 1) {
            new_individus[i] = joueurs->individus[classement1[i]];
        }
    }
    else {
        for (int i = 0; i != joueurs->capacite; i += 1) {
            new_individus[i] = joueurs->individus[classement2[i]];
        }
    }
    free(classement1);
    free(classement2);
    free(coeffs_binomiaux);
    free(joueurs->individus);
    joueurs->individus = new_individus;
}

void tournoi2(pop2* mutateurs, int nb_individus) {
    for (int i = 0; i != mutateurs->capacite; i += 1) {
        mut* mutateur = mutateurs->individus[i];
        mutateur->score = INT_MIN;
        for (int j = 0; j != mutateur->nb_enfants; j += 1) {
            mutateur->score = max(mutateur->score, mutateur->enfants[j]->score);
        }
    }
    mut** new_individus = malloc(sizeof(mut*) * mutateurs->capacite);
    int c = 0;
    for (int i = nb_individus; i >= -nb_individus; i -= 1) {
        for (int j = 0; j != mutateurs->capacite; j += 1) {
            if (mutateurs->individus[j]->score == i) {
                new_individus[c] = mutateurs->individus[j];
                c += 1;
            }
        }
    }
    free(mutateurs->individus);
    mutateurs->individus = new_individus;
}

ind* creer_ind(mut* mutateur) {
    double* data = malloc(sizeof(double) * mutateur->nb_params);
    for (int i = 0; i != mutateur->nb_params; i += 1) {
        data[i] = rand_normal(1) * mutateur->data[i];
    }
    ind* joueur = malloc(sizeof(ind));
    joueur->data = data;
    joueur->score = 0;
    joueur->pere = mutateur;
    joueur->nb_params = mutateur->nb_params;
    mutateur->nb_enfants_en_vie += 1;
    return joueur;
}

mut* creer_mut(int nb_enfants, int nb_params) {
    double* data = malloc(sizeof(double) * nb_params);
    ind** enfants = malloc(sizeof(ind*) * nb_enfants);
    for (int i = 0; i != nb_params; i += 1) {
        data[i] = exp(rand_normal(1));
    }
    mut* mutateur = malloc(sizeof(mut));
    mutateur->data = data;
    mutateur->nb_enfants = nb_enfants;
    mutateur->enfants = enfants;
    mutateur->score = INT_MIN;
    mutateur->nb_enfants_en_vie = 0;
    mutateur->nb_params = nb_params;
    return mutateur;
}

void muter_ind(ind* joueur, mut* mutateur) {
    for (int i = 0; i != joueur->nb_params; i += 1) {
        joueur->data[i] += rand_normal(mutateur->data[i]);
    }
    joueur->pere = mutateur;
    mutateur->nb_enfants_en_vie += 1;
}

void muter_mut(mut* mutateur) {
    for (int i = 0; i != mutateur->nb_params; i += 1) {
        mutateur->data[i] *= exp(rand_normal(1));
    }
}

ind* copy_ind(ind* joueur) {
    double* data = malloc(sizeof(double) * joueur->nb_params);
    for (int i = 0; i != joueur->nb_params; i += 1) {
        data[i] = joueur->data[i];
    }
    ind* new_joueur = malloc(sizeof(ind));
    new_joueur->data = data;
    new_joueur->score = 0;
    new_joueur->pere = NULL;
    new_joueur->nb_params = joueur->nb_params;
    return new_joueur;
}

mut* copy_mut(mut* mutateur) {
    double* data = malloc(sizeof(double) * mutateur->nb_params);
    ind** enfants = malloc(sizeof(ind*) * mutateur->nb_enfants);
    for (int i = 0; i != mutateur->nb_params; i += 1) {
        data[i] = mutateur->data[i];
    }
    mut* new_mutateur = malloc(sizeof(mut));
    new_mutateur->data = data;
    new_mutateur->nb_enfants = mutateur->nb_enfants;
    new_mutateur->enfants = enfants;
    new_mutateur->score = INT_MIN;
    new_mutateur->nb_enfants_en_vie = 0;
    new_mutateur->nb_params = mutateur->nb_params;
    return new_mutateur;
}

pop2* creer_pop2(int nb_mutateurs, int nb_enfants, int nb_params) {
    int capacite = 1 << nb_mutateurs;
    mut** individus = malloc(sizeof(mut*) * capacite);
    for (int i = 0; i != capacite; i += 1) {
        individus[i] = creer_mut(nb_enfants, nb_params);
    }
    pop2* mutateurs = malloc(sizeof(pop2));
    mutateurs->individus = individus;
    mutateurs->capacite = capacite;
    mutateurs->log2_capacite = nb_mutateurs;
    return mutateurs;
}

pop* creer_pop(int nb_joueurs, int nb_mutateurs, int nb_params) {
    int capacite = 1 << nb_joueurs;
    int capacite2 = 1 << nb_mutateurs;
    int nb_enfants = 1 << (nb_joueurs - nb_mutateurs - 1);
    ind** individus = malloc(sizeof(ind*) * capacite);
    pop2* mutateurs = creer_pop2(nb_mutateurs, nb_enfants, nb_params);
    for (int i = 0; i != capacite2; i += 1) {
        for (int j = 0; j != nb_enfants; j += 1) {
            individus[nb_enfants * i + j] = creer_ind(mutateurs->individus[i]);
            mutateurs->individus[i]->enfants[j] = individus[nb_enfants * i + j];
        }  
    }
    for (int i = 0; i != capacite2; i += 1) {
        for (int j = 0; j != nb_enfants; j += 1) {
            individus[nb_enfants * i + j + capacite / 2] = creer_ind(mutateurs->individus[i]);
        }  
    }
    pop* joueurs = malloc(sizeof(pop));
    joueurs->individus = individus;
    joueurs->capacite = capacite;
    joueurs->log2_capacite = nb_joueurs;
    joueurs->mutateurs = mutateurs;
    joueurs->nb_vic = 0;
    return joueurs;
}

pops* creer_pops(int nb_params) {
    int nb_joueurs = 8;
    int nb_mutateurs = 4;
    pop** individus = malloc(9 * sizeof(pop*));
    for (int i = 0; i != 3; i += 1) {
        for (int j = 0; j != 3; j += 1) {
            individus[3 * i + j] = creer_pop(nb_joueurs - 1 + i, nb_mutateurs - 2 + i + j, nb_params);
        }
    }
    int* victoires = malloc(3 * sizeof(int));
    for (int i = 0; i != 3; i += 1) {
        victoires[i] = 0;
    }
    pops* populations = malloc(sizeof(pops));
    populations->individus = individus;
    populations->plus_faible = 4;
    populations->victoires = victoires;
    populations->nb_params = nb_params;
    return populations;
}

void lib_mem_mut(mut* mutateur) {
    free(mutateur->data);
    free(mutateur->enfants);
    free(mutateur);
}

void diminuer(mut* mutateur) {
    mutateur->nb_enfants_en_vie -= 1;
    if (mutateur->nb_enfants_en_vie == -1) {
        lib_mem_mut(mutateur);
    }
}

void lib_mem_ind(ind* joueur) {
    free(joueur->data);
    if (joueur->pere != NULL) {
        diminuer(joueur->pere);
    }  
    free(joueur);
}

void lib_mem_pop2(pop2* mutateurs) {
    if (mutateurs != NULL) {
        for (int i = 0; i != mutateurs->capacite; i += 1) {
            diminuer(mutateurs->individus[i]);
        }
        free(mutateurs->individus);
        free(mutateurs);
    }
}

void lib_mem_pop(pop* joueurs) {
    for (int i = 0; i != joueurs->capacite; i += 1) {
        lib_mem_ind(joueurs->individus[i]);
    }
    lib_mem_pop2(joueurs->mutateurs);
    free(joueurs->individus);
    free(joueurs);
}

void lib_mem_pops(pops* populations) {
    for (int i = 0; i != 9; i += 1) {
        lib_mem_pop(populations->individus[i]);
    }
    free(populations->individus);
    free(populations->victoires);
    free(populations);
}

void selection(pop* joueurs) {
    for (int i = joueurs->capacite / 2; i != joueurs->capacite; i += 1) {
        lib_mem_ind(joueurs->individus[i]);
        joueurs->individus[i] = NULL;
    }
}

void selection2(pop2* mutateurs) {
    for (int i = mutateurs->capacite / 2; i != mutateurs->capacite; i += 1) {
        diminuer(mutateurs->individus[i]);
        mutateurs->individus[i] = NULL;
    }
}

void reproduction(pop* joueurs) {
    int nb = joueurs->mutateurs->capacite;
    int nb2 = joueurs->capacite / 2;
    int nb_enfants = nb2 / nb;
    for (int i = 0; i != nb; i += 1) {
        for (int j = 0; j != nb_enfants; j += 1) {
            joueurs->individus[i + nb * j + nb2] = copy_ind(joueurs->individus[i + nb * j]);
            muter_ind(joueurs->individus[i + nb * j + nb2], joueurs->mutateurs->individus[i]);
            joueurs->mutateurs->individus[i]->enfants[j] = joueurs->individus[i + nb * j + nb2];
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

void evolution(pop* joueurs, bool (*f) (double*, double*)) {
    tournoi(joueurs, f);
    tournoi2(joueurs->mutateurs, joueurs->log2_capacite);
    selection2(joueurs->mutateurs);
    selection(joueurs);
    reproduction2(joueurs->mutateurs);
    reproduction(joueurs);
}

void diviser_pop2(pop2* mutateurs) {
    int nb_enfants = mutateurs->individus[0]->nb_enfants * 2;
    for (int i = mutateurs->capacite / 4; i != mutateurs->capacite; i += 1) {
        diminuer(mutateurs->individus[i]);
        mutateurs->individus[i] = NULL;
    }
    mutateurs->capacite /= 2;
    mutateurs->log2_capacite -= 1;
    mut** new_individus = malloc(sizeof(mut*) * mutateurs->capacite);
    for (int i = 0; i != mutateurs->capacite / 2; i += 1) {
        ind** new_enfants = malloc(sizeof(ind*) * nb_enfants);
        free(mutateurs->individus[i]->enfants);
        new_individus[i] = mutateurs->individus[i];
        new_individus[i]->enfants = new_enfants;
        new_individus[i]->nb_enfants = nb_enfants;
    }
    free(mutateurs->individus);
    mutateurs->individus = new_individus;
    reproduction2(mutateurs);
}

void multiplier_pop2(pop2* mutateurs) {
    int nb_enfants = mutateurs->individus[0]->nb_enfants / 2;
    mutateurs->capacite *= 2;
    mutateurs->log2_capacite += 1;
    mut** new_individus = malloc(sizeof(mut*) * mutateurs->capacite);
    for (int i = 0; i != mutateurs->capacite / 2; i += 1) {
        ind** new_enfants = malloc(sizeof(ind*) * nb_enfants);
        free(mutateurs->individus[i]->enfants);
        new_individus[i] = mutateurs->individus[i];
        new_individus[i]->enfants = new_enfants;
        new_individus[i]->nb_enfants = nb_enfants;
    }
    free(mutateurs->individus);
    mutateurs->individus = new_individus;
    reproduction2(mutateurs);
}

void diviser_pop(pop* joueurs) {
    for (int i = joueurs->capacite / 4; i != joueurs->capacite; i += 1) {
        lib_mem_ind(joueurs->individus[i]);
        joueurs->individus[i] = NULL;
    }
    joueurs->capacite /= 2;
    ind** new_individus1 = malloc(sizeof(ind*) * joueurs->capacite);
    for (int i = 0; i != joueurs->capacite / 2; i += 1) {
        new_individus1[i] = joueurs->individus[i];
    }
    free(joueurs->individus);
    joueurs->individus = new_individus1;
    joueurs->log2_capacite -= 1;
    for (int i = joueurs->mutateurs->capacite / 4; i != joueurs->mutateurs->capacite; i += 1) {
        diminuer(joueurs->mutateurs->individus[i]);
        joueurs->mutateurs->individus[i] = NULL;
    }
    joueurs->mutateurs->capacite /= 2;
    joueurs->mutateurs->log2_capacite -= 1;
    mut** new_individus2 = malloc(sizeof(mut*) * joueurs->mutateurs->capacite);
    for (int i = 0; i != joueurs->mutateurs->capacite / 2; i += 1) {
        new_individus2[i] = joueurs->mutateurs->individus[i];
    }
    free(joueurs->mutateurs->individus);
    joueurs->mutateurs->individus = new_individus2;
    reproduction2(joueurs->mutateurs);
    reproduction(joueurs);
}

void multiplier_pop(pop* joueurs) {
    joueurs->capacite *= 2;
    ind** new_individus1 = malloc(sizeof(ind*) * joueurs->capacite);
    for (int i = 0; i != joueurs->capacite / 2; i += 1) {
        new_individus1[i] = joueurs->individus[i];
    }
    free(joueurs->individus);
    joueurs->individus = new_individus1;
    joueurs->log2_capacite += 1;
    joueurs->mutateurs->capacite *= 2;
    joueurs->mutateurs->log2_capacite += 1;
    mut** new_individus2 = malloc(sizeof(mut*) * joueurs->mutateurs->capacite);
    for (int i = 0; i != joueurs->mutateurs->capacite / 2; i += 1) {
        new_individus2[i] = joueurs->mutateurs->individus[i];
    }
    free(joueurs->mutateurs->individus);
    joueurs->mutateurs->individus = new_individus2;
    reproduction2(joueurs->mutateurs);
    reproduction(joueurs);
}

void evolutions(pops* populations, bool (*f)(double*, double*)) {
    for (int i = 0; i != 3; i += 1) {
        for (int j = 0; j != 3; j += 1) {          
            for (int k = 0; k != (1 << (2 - i)) - 1; k += 1) {
                evolution(populations->individus[3 * i + j], f);
            }
            tournoi(populations->individus[3 * i + j], f);
            tournoi2(populations->individus[3 * i + j]->mutateurs, populations->individus[3 * i + j]->log2_capacite);
        }
    }
    ind** individus = malloc(sizeof(ind*) * 8);
    int j = 0;
    for (int i = 0; i != 9; i += 1) {
        if (i != populations->plus_faible) {
            individus[j] = copy_ind(populations->individus[i]->individus[0]);
            j += 1;
        }
    }
    ind** mem = malloc(sizeof(ind*) * 9);
    for (int i = 0; i != 9; i += 1) {
        if (i != populations->plus_faible) {
            mem[i] = individus[i - ((i < populations->plus_faible) ? 0 : 1)];
        }
        else {
            mem[i] = NULL;
        }
    }
    pop* best_players = malloc(sizeof(pop));
    best_players->individus = individus;
    best_players->mutateurs = NULL;
    best_players->capacite = 8;
    best_players->log2_capacite = 3;
    tournoi(best_players, f);
    int premier = 0;
    int dernier = 0;
    int premier1 = 0;
    int premier2 = 3;
    int premier3 = 6;
    int i = 0;
    individus = best_players->individus;
    while (individus[i] != mem[0] && individus[i] != mem[1] && individus[i] != mem[2]) {
        i += 1;
    }
    while (individus[i] != mem[premier1]) {
        premier1 += 1;
    }
    i = 0;
    while (individus[i] != mem[3] && individus[i] != mem[4] && individus[i] != mem[5]) {
        i += 1;
    }
    while (individus[i] != mem[premier2]) {
        premier2 += 1;
    }
    i = 0;
    while (individus[i] != mem[6] && individus[i] != mem[7] && individus[i] != mem[8]) {
        i += 1;
    }
    while (individus[i] != mem[premier3]) {
        premier3 += 1;
    }
    while (mem[premier] != individus[0]) {
        premier += 1;
    }
    while (mem[dernier] != individus[7]) {
        dernier += 1;
    }
    for (int i = 0; i != 9; i += 1) {
        if (populations->individus[i]->individus[0] != populations->individus[premier]->individus[0]) {
            lib_mem_ind(populations->individus[i]->individus[populations->individus[i]->capacite / 4 - 1]);
            populations->individus[i]->individus[populations->individus[i]->capacite / 4 - 1] = copy_ind(populations->individus[premier]->individus[0]);
            populations->individus[i]->individus[populations->individus[i]->capacite / 4 - 1]->pere = populations->individus[premier]->individus[0]->pere;
            populations->individus[premier]->individus[0]->pere->nb_enfants_en_vie += 1;
        }
        if (populations->individus[i]->individus[0]->pere != populations->individus[premier]->individus[0]->pere && populations->individus[i]->mutateurs->log2_capacite > 2) {
            mut* copie = copy_mut(populations->individus[i]->mutateurs->individus[populations->individus[i]->mutateurs->capacite / 4 - 1]);
            diminuer(populations->individus[i]->mutateurs->individus[populations->individus[i]->mutateurs->capacite / 4 - 1]);
            for (int j = 0; j != populations->nb_params; j += 1) {
                copie->data[j] = populations->individus[premier]->individus[0]->pere->data[j];
            }
            populations->individus[i]->mutateurs->individus[populations->individus[i]->mutateurs->capacite / 4 - 1] = copie;
        }
    }
    populations->plus_faible = dernier;
    bool pass = true;
    for (int i = 0; i != 9; i += 3) {
        if (premier == i || premier == i + 1 || premier == i + 2) {
            populations->victoires[i / 3] += 1;
            if (populations->victoires[i / 3] >= 3) {
                if (i / 3 == 0 && populations->individus[0]->mutateurs->log2_capacite > 2 && populations->individus[3]->mutateurs->log2_capacite > 2 && populations->individus[6]->mutateurs->log2_capacite > 2) {
                    for (int j = 0; j != 9; j += 1) {
                        diviser_pop(populations->individus[j]);
                    }
                    pass = false;
                }
                if (i / 3 == 2) {
                    for (int j = 0; j != 9; j += 1) {
                        multiplier_pop(populations->individus[j]);
                    }
                    pass = false;
                } 
                populations->victoires[i / 3] = 0; 
            }
        }
        else {
            populations->victoires[i / 3] = 0;
        }
    }
    bool pass1 = true;
    bool pass2 = true;
    bool pass3 = true;
    if (pass) {
        for (int i = 0; i != 3; i += 1) {
            if (premier1 == i) {
                populations->individus[i]->nb_vic += 1;
                if (populations->individus[i]->nb_vic >= 5) {
                    if (i == 0 && populations->individus[0]->mutateurs->log2_capacite > 2) {
                        for (int j = 0; j != 3; j += 1) {
                            diviser_pop2(populations->individus[j]->mutateurs);
                        }
                        pass1 = false;
                    }
                    if (i == 2 && populations->individus[2]->mutateurs->log2_capacite + 1 < populations->individus[2]->log2_capacite) {
                        for (int j = 0; j != 3; j += 1) {
                            multiplier_pop2(populations->individus[j]->mutateurs);
                        }
                        pass1 = false;
                    }
                    populations->individus[i]->nb_vic = 0;
                }
            }
            else {
                populations->individus[i]->nb_vic = 0;
            }
        }
        for (int i = 3; i != 6; i += 1) {
            if (premier2 == i) {
                populations->individus[i]->nb_vic += 1;
                if (populations->individus[i]->nb_vic >= 5) {
                    if (i == 3 && populations->individus[3]->mutateurs->log2_capacite > 2) {
                        for (int j = 3; j != 6; j += 1) {
                            diviser_pop2(populations->individus[j]->mutateurs);
                        }
                        pass2 = false;
                    }
                    if (i == 5 && populations->individus[5]->mutateurs->log2_capacite + 1 < populations->individus[5]->log2_capacite) {
                        for (int j = 3; j != 6; j += 1) {
                            multiplier_pop2(populations->individus[j]->mutateurs);
                        }
                        pass2 = false;
                    }
                    populations->individus[i]->nb_vic = 0;
                }
            }
            else {
                populations->individus[i]->nb_vic = 0;
            }
        }
        for (int i = 6; i != 9; i += 1) {
            if (premier3 == i) {
                populations->individus[i]->nb_vic += 1;
                if (populations->individus[i]->nb_vic >= 5) {
                    if (i == 6 && populations->individus[6]->mutateurs->log2_capacite > 2) {
                        for (int j = 6; j != 9; j += 1) {
                            diviser_pop2(populations->individus[j]->mutateurs);
                        }
                        pass3 = false;
                    }
                    if (i == 8 && populations->individus[8]->mutateurs->log2_capacite + 1 < populations->individus[8]->log2_capacite) {
                        for (int j = 6; j != 9; j += 1) {
                            multiplier_pop2(populations->individus[j]->mutateurs);
                        }
                        pass3 = false;
                    }
                    populations->individus[i]->nb_vic = 0;
                }
            }
            else {
                populations->individus[i]->nb_vic = 0;
            }
        }
    }
    if (pass) {
        if (pass1) {
            for (int i = 0; i != 3; i += 1) {
                selection2(populations->individus[i]->mutateurs);
                selection(populations->individus[i]);
                reproduction2(populations->individus[i]->mutateurs);
                reproduction(populations->individus[i]);
            }
        }
        else {
            for (int i = 0; i != 3; i += 1) {
                selection(populations->individus[i]);
                reproduction(populations->individus[i]);
            }
        }
        if (pass2) {
            for (int i = 3; i != 6; i += 1) {
                selection2(populations->individus[i]->mutateurs);
                selection(populations->individus[i]);
                reproduction2(populations->individus[i]->mutateurs);
                reproduction(populations->individus[i]);
            }
        }
        else {
            for (int i = 3; i != 6; i += 1) {
                selection(populations->individus[i]);
                reproduction(populations->individus[i]);
            }
        }
        if (pass3) {
            for (int i = 6; i != 9; i += 1) {
                selection2(populations->individus[i]->mutateurs);
                selection(populations->individus[i]);
                reproduction2(populations->individus[i]->mutateurs);
                reproduction(populations->individus[i]);
            }
        }
        else {
            for (int i = 6; i != 9; i += 1) {
                selection(populations->individus[i]);
                reproduction(populations->individus[i]);
            }
        }
    }
    free(mem);
    lib_mem_pop(best_players);
}

void maj_parametres(pops* populations, bool (*f)(double*, double*), double* parametres) {
    int i = 0;
    for (int j = 1; j != 9; j += 1) {
        if (f(populations->individus[j]->individus[0]->data, populations->individus[i]->individus[0]->data)) {
            i = j;
        }
    }
    for (int j = 0; j != populations->nb_params; j += 1) {
        parametres[j] = populations->individus[i]->individus[0]->data[j];
    }
}

pops* populations = NULL;

void initialiser(int nb_parametres) {
    populations = creer_pops(nb_parametres);
}
