#include "eval.h"
#include <math.h>
#include <cassert>
#include <armadillo>

using namespace std;
double e_Au = 2.951;
double r_Au =5.29; 
// double r_Au =2.0; 
// double e_Au = 2;

void ConvertAtomstoCoor(const vector<Atom> &Atoms, arma::mat &Coor){
    for(size_t k = 0; k < Atoms.size(); k++)
        for (size_t j = 0; j < 3; j++)
            Coor(j,k) = Atoms[k].r[j];
}

void ConvertCoortoAtoms(vector<Atom> &Atoms, const arma::mat &Coor){
    for(size_t k = 0; k < Atoms.size(); k++)
        for (size_t j = 0; j < 3; j++)
            Atoms[k].r[j] = Coor(j,k);
}

double E_LJ(const std::vector<Atom> &Atoms)
{
    double E = 0.0;
    for (size_t k = 0; k < Atoms.size(); k++){
        Atom atom_k = Atoms[k];
        for (size_t j = k + 1; j < Atoms.size(); j++)
        {
            Atom atom_j = Atoms[j];
            double R2 = (atom_j.r[0] - atom_k.r[0]) * (atom_j.r[0] - atom_k.r[0])
                 + (atom_j.r[1] - atom_k.r[1]) * (atom_j.r[1] - atom_k.r[1]) 
                 + (atom_j.r[2] - atom_k.r[2]) * (atom_j.r[2] - atom_k.r[2]);
            double r2overR2 =  r_Au * r_Au / R2;
            double E_temp = pow(r2overR2, 6) - 2.0 * pow(r2overR2, 3);
            E += E_temp * e_Au;
            // cout << E_temp * e_Au << endl;
        }
    }
    return E;
}

double E_LJ(const arma::mat &coords) {
    int n_atoms = coords.n_cols;
    std::vector<Atom> Atoms(n_atoms);
    ConvertCoortoAtoms(Atoms, coords);
    // printf("current energy: %.3e\n", E_LJ(Atoms));
    return E_LJ(Atoms);
}


void F_LJ_forward(arma::mat &force, const std::vector<Atom> &Atoms, double stepsize){
    vector<Atom> Atoms_forward = Atoms;
    assert(force.n_rows == 3 && force.n_cols == Atoms.size());
    double E_unmove = E_LJ(Atoms);
    for (size_t k = 0; k < Atoms.size(); k++) {
        for (size_t j = 0; j < 3; j++) {
            Atoms_forward = Atoms;
            Atoms_forward[k].r[j] += stepsize;
            double E_forward = E_LJ(Atoms_forward);
            force(j, k) = -(E_forward - E_unmove)/stepsize;
        }
    }
}

void F_LJ_central(arma::mat &force, const vector<Atom> &Atoms, double stepsize)
{
    vector<Atom> Atoms_forward = Atoms;
    vector<Atom> Atoms_backward = Atoms;
    assert(force.n_rows == 3 && force.n_cols == Atoms.size());
    for (size_t k = 0; k < Atoms.size(); k++){
        for (size_t j = 0; j < 3; j++){
            Atoms_forward = Atoms;
            Atoms_backward = Atoms;
            Atoms_forward[k].r[j] = Atoms_forward[k].r[j] + stepsize;
            Atoms_backward[k].r[j] = Atoms_backward[k].r[j] - stepsize;
            // cout << Atoms_forward[k].r[j] << endl;
            // cout << Atoms_backward[k].r[j] << endl;
            double E_forward = E_LJ(Atoms_forward);
            double E_backward = E_LJ(Atoms_backward);
            // cout << E_forward << endl;
            // cout << E_backward << endl;
            force(j,k) = -(E_forward - E_backward) /2.0 /stepsize;
        }
    }
}

void F_LJ_analytic(arma::mat &force, const std::vector<Atom> &Atoms){
    for (int k = 0; k < Atoms.size(); k++) {
        for (int j = 0; j < 3; j++) {
            double f = 0.0;
            for (int i = 0; i < Atoms.size(); i++) {
                if (i != k) {
                    double Rik = sqrt(pow(Atoms[k].r[0]-Atoms[i].r[0],2)+
                    pow(Atoms[k].r[1]-Atoms[i].r[1],2)+pow(Atoms[k].r[2]-Atoms[i].r[2],2));
                    f += 12*e_Au*(Atoms[i].r[j]-Atoms[k].r[j])*(pow(r_Au,6)/pow(Rik,8)-pow(r_Au,12)/pow(Rik,14));
                }
            }
            force(j,k) = f;
        }
    }
}

void Steepest_descend( vector<Atom> &opt_Atoms, const vector<Atom> &Atoms, double fdstepsize, double thresh){

    double search_stepsize = 0.3;
    double E_old = E_LJ(Atoms), E_new;
    vector<Atom> new_Atoms = Atoms;
    arma::mat new_point(3, Atoms.size());
    arma::mat old_point(3, Atoms.size());
    arma::mat Force(3, Atoms.size());
    arma::mat fForce(3, Atoms.size());
    arma::mat aForce(3, Atoms.size());
    ConvertAtomstoCoor(Atoms, old_point);

    // evaluate the Force at starting point
    F_LJ_forward(fForce, Atoms, fdstepsize);
    fForce.print("Forward fd Force");
    cout << "Forward fd Force norm: "<< arma::norm(fForce, "fro") << endl;
    F_LJ_central(Force, Atoms, fdstepsize);
    Force.print("Central fd Force");
    cout << "Central fd Force norm: "<< arma::norm(Force, "fro") << endl;
    F_LJ_analytic(aForce, Atoms);
    aForce.print("Analytic Force");
    cout << "Analytic Force norm: "<< arma::norm(aForce, "fro") << endl;
    cout << "Energy: "<< E_LJ(Atoms) << endl;
    int count = 0;
    cout << endl << "Start opt with central finite difference force evaluation" << endl << endl;
    while (arma::norm(Force, "fro") > thresh && count < 1e2){
        // calculate new point position
        new_point = old_point + search_stepsize * Force / arma::norm(Force, 2);
        ConvertCoortoAtoms(new_Atoms, new_point);
        E_new = E_LJ(new_Atoms);
        cout << "count "<< count  << " Energy: "<< E_new << endl;
        // the step makes function evaluation lower - it is a good step. what do you do?
        if (E_new < E_old){
            old_point = new_point;
            E_old = E_new;
            // refresh Force
            F_LJ_central(Force, new_Atoms, fdstepsize);
            search_stepsize *= 1.2;
            new_point.print("new_point");
            cout << "Force norm: "<< arma::norm(Force, "fro") << endl;
        }
        else // the step makes function evaluation higher - it is a bad step. what do you do?
            search_stepsize /= 2;
        count += 1;

    }

    cout << "Total iterations: " << count << endl;
    printf("Final energy: %.3e\n", E_LJ(old_point));
    ConvertCoortoAtoms(opt_Atoms, old_point);
}

void Steepest_descend_line_search( std::vector<Atom> &opt_Atoms, const std::vector<Atom> &Atoms,
double fdstepsize, double thresh) {
    double search_stepsize = 0.3;
    double golden_ratio = 0.38197;
    // first bracket A<B<D (A is the inital point), such that E(B)<E(A),E(B)<E(D)
    // during line search, we keep the order A<B<C<D
    arma::mat A_point(3, Atoms.size());
    ConvertAtomstoCoor(Atoms, A_point);
    arma::mat B_point(3, Atoms.size());
    arma::mat C_point(3, Atoms.size());
    arma::mat D_point(3, Atoms.size());
    // length between points
    double AB, BD, AD;
    double golden_tole = 1e-7; // tolerance to stop golden search
    // forces
    arma::mat Force(3, Atoms.size());
    arma::mat unit_Force(3, Atoms.size());
    arma::mat fForce(3, Atoms.size());
    arma::mat aForce(3, Atoms.size());

    // evaluate the Force at starting point
    F_LJ_forward(fForce, Atoms, fdstepsize);
    fForce.print("Forward fd Force");
    cout << "Forward fd Force norm: "<< arma::norm(fForce, "fro") << endl;
    F_LJ_central(Force, Atoms, fdstepsize);
    Force.print("Central fd Force");
    cout << "Central fd Force norm: "<< arma::norm(Force, "fro") << endl;
    F_LJ_analytic(aForce, Atoms);
    aForce.print("Analytic Force");
    cout << "Analytic Force norm: "<< arma::norm(aForce, "fro") << endl;
    cout << "Energy: "<< E_LJ(Atoms) << endl;
    int count = 0;
    cout << endl << "Start opt with central finite difference force evaluation" << endl << endl;
    while (arma::norm(Force, "fro") > thresh && count < 5e2){
        unit_Force = Force / arma::norm(Force, 2);
        // get B_point
        int bracket_count = 0;
        B_point = A_point + search_stepsize * unit_Force;
        while (E_LJ(B_point) > E_LJ(A_point)) {
            bracket_count++;
            search_stepsize /= 2;
            B_point = A_point + search_stepsize * unit_Force;
            printf("Bracketing for line search, finding B point with stepsize %.2f\n", search_stepsize);
            if (bracket_count > 100) {
                printf("Finding B point failed, please adjust search_stepsize\n");
                goto step_forward;
            }
        }
        cout << "Bracketing for line search, B point found" << endl;
        // get D_point
        bracket_count = 0;
        D_point = B_point + search_stepsize * unit_Force;
        while (E_LJ(D_point) < E_LJ(B_point)) {
            bracket_count++;
            search_stepsize *= 1.2;
            D_point = B_point + search_stepsize * unit_Force;
            printf("Bracketing for line search, finding D point with stepsize %.2f\n", search_stepsize);
            if (bracket_count > 100) {
                printf("Finding D point failed, please adjust search_stepsize\n");
                goto step_forward;
            }
        }
        cout << "Bracketing for line search, D point found" << endl;
        // if we fail to find B or D, it is possible there's no local minima along the
        // negative gradient direction, in this case we move one step towards the negative
        // gradient to decrease energy, then use line search in the future steps
        step_forward:
        if(bracket_count > 100) {
            printf("Bracketing failed, use normal steepest descend instead\n");
            search_stepsize = 0.3;
            arma::mat new_point = A_point + search_stepsize*unit_Force;
            while (E_LJ(new_point) > E_LJ(A_point)) {
                search_stepsize /= 2;
            }
            A_point = new_point;
            ConvertCoortoAtoms(opt_Atoms, A_point);
            printf("current energy: %.3e\n", E_LJ(A_point));
            A_point.print("new_point");
            goto update_force;
        }
        // Golden search
        AB = arma::norm(B_point-A_point, "fro");
        BD = arma::norm(D_point-B_point, "fro");
        AD = arma::norm(D_point-A_point, "fro");
        if (AB < BD) {
            C_point = D_point + golden_ratio*(A_point-D_point);
        } else {
            C_point = B_point;
            B_point = A_point + golden_ratio*(D_point-A_point);
        }
        while (AD > golden_tole) {
            if (E_LJ(B_point) > E_LJ(C_point)) {
                A_point = B_point;
                B_point = C_point;
            } else {
                D_point = C_point;
            }
            AB = arma::norm(B_point-A_point, "fro");
            BD = arma::norm(D_point-B_point, "fro");
            AD = arma::norm(D_point-A_point, "fro");
            if (AB < BD) {
                C_point = D_point + golden_ratio*(A_point-D_point);
            } else {
                C_point = B_point;
                B_point = A_point + golden_ratio*(D_point-A_point);
            }
        }
        if (E_LJ(B_point) > E_LJ(C_point)) {
            ConvertCoortoAtoms(opt_Atoms, C_point);
            printf("current energy: %.3e\n", E_LJ(C_point));
            C_point.print("new_point");
        } else {
            ConvertCoortoAtoms(opt_Atoms, B_point);
            printf("current energy: %.3e\n", E_LJ(B_point));
            B_point.print("new_point");
        }
        // revaluate force in the new position
        update_force:
        F_LJ_analytic(Force, opt_Atoms);
        // F_LJ_central(Force, opt_Atoms, fdstepsize);
        cout << "Force norm: "<< arma::norm(Force, "fro") << endl;
        count += 1;
    }

    cout << "Total iterations: " << count << endl;
    printf("Final energy: %.3e\n", E_LJ(opt_Atoms));
}
