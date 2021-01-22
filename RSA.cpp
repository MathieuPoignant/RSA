#include <chrono>
#include <iostream>
#include <random>
using namespace std;

class GrandNombre {
 public:
  GrandNombre();
  ~GrandNombre();
  GrandNombre(int size);
  GrandNombre(int size, const unsigned long long int *value);
  GrandNombre resize(int new_size) const;
  void Print() const;
  void operator=(const GrandNombre &copied);
  GrandNombre operator+(const GrandNombre &added) const;
  GrandNombre operator-(const GrandNombre &substracted) const;
  GrandNombre operator*(const GrandNombre &multiplied) const;
  bool operator<=(const GrandNombre &compared) const;
  GrandNombre operator<<(int decalage) const;
  GrandNombre operator>>(int decalage) const;
  GrandNombre operator%(const GrandNombre &modulo) const;
  GrandNombre subMod(const GrandNombre &substracted,
                     const GrandNombre &modulo) const;
  GrandNombre Montgomery(const GrandNombre &multiplied, const GrandNombre &R,
                         const GrandNombre &N, const GrandNombre &V,
                         const GrandNombre &invV) const;
  GrandNombre Phi(const GrandNombre &R, const GrandNombre &N) const;
  GrandNombre Pow(const GrandNombre &power, const GrandNombre &R,
                  const GrandNombre &N, const GrandNombre &V,
                  const GrandNombre &invV) const;
  GrandNombre friend RSA(const GrandNombre &M, const GrandNombre &E,
                         const GrandNombre &D, const GrandNombre &R,
                         const GrandNombre &N, const GrandNombre &V,
                         const GrandNombre &invV);

 private:
  int tab_size;
  unsigned long long int *tab;
};

GrandNombre::GrandNombre() {
  tab_size = 0;
  unsigned long long int *newtab = new unsigned long long int[tab_size];
  tab = newtab;
}
GrandNombre::~GrandNombre() { delete[] tab; }

GrandNombre::GrandNombre(int size) {
  tab_size = size;
  unsigned long long int *newtab = new unsigned long long int[tab_size];
  for (int i = 0; i < tab_size; i++) {
    newtab[i] = 0;
  }
  tab = newtab;
}

GrandNombre::GrandNombre(int size, const unsigned long long int *value) {
  tab_size = size;
  unsigned long long int *newtab = new unsigned long long int[tab_size];
  for (int i = 0; i < tab_size; i++) {
    newtab[i] = value[i];
  }
  tab = newtab;
}

GrandNombre GrandNombre::resize(int new_size) const {
  GrandNombre newGrandNombre = GrandNombre(new_size);
  int decalage = this->tab_size - new_size;
  if (decalage >= 0) {  // enlève les poids forts (0)
    for (int i = 0; i < new_size; i++) {
      newGrandNombre.tab[i] = this->tab[i + decalage];
    }
  } else {  // rajoute des 0 en poids forts
    for (int i = 0; i < 0 - decalage; i++) {
      newGrandNombre.tab[i] = 0;
    }
    for (int i = 0; i < this->tab_size; i++) {
      newGrandNombre.tab[i - decalage] = this->tab[i];
    }
  }
  return newGrandNombre;
}

void GrandNombre::Print() const {
  cout << "[";
  for (int i = 0; i < this->tab_size; i++) {
    cout << tab[i] << ", ";
  }
  cout << "]" << endl;
}

void GrandNombre::operator=(const GrandNombre &copied) {
  // recopie les poids faibles (souvent des 0 en poids fort)
  for (int i = this->tab_size; i >= 0; i--) {
    this->tab[this->tab_size - i] = copied.tab[copied.tab_size - i];
  }
}

GrandNombre GrandNombre::operator+(const GrandNombre &added) const {
  int size = max(this->tab_size, added.tab_size);
  GrandNombre newGrandNombre = GrandNombre(size);
  unsigned long long int temp;
  unsigned long long int carry = 0;
  int decalage = this->tab_size - added.tab_size;
  if (decalage >= 0) {  // this plus grand
    for (int i = added.tab_size - 1; i >= 0; i--) {
      temp = this->tab[i + decalage] + added.tab[i] + carry;
      carry = temp / 4294967296;
      newGrandNombre.tab[i + decalage] = temp % 4294967296;
    }
    // recopie les poids forts de this en propageant la retenue
    for (int i = decalage - 1; i >= 0; i--) {
      temp = this->tab[i] + carry;
      carry = temp / 4294967296;
      newGrandNombre.tab[i] = temp % 4294967296;
    }
  } else {  // this plus petit
    for (int i = this->tab_size - 1; i >= 0; i--) {
      temp = this->tab[i] + added.tab[i - decalage] + carry;
      carry = temp / 4294967296;
      newGrandNombre.tab[i - decalage] = temp % 4294967296;
    }
    // recopie les poids forts de added en propageant la retenue
    for (int i = -decalage - 1; i >= 0; i--) {
      temp = added.tab[i] + carry;
      carry = temp / 4294967296;
      newGrandNombre.tab[i] = temp % 4294967296;
    }
  }
  return newGrandNombre;
}

GrandNombre GrandNombre::operator-(const GrandNombre &substracted) const {
  int size = max(this->tab_size, substracted.tab_size);
  GrandNombre newGrandNombre = GrandNombre(size);
  unsigned long long int temp;
  unsigned long long int carry = 0;
  int decalage = this->tab_size - substracted.tab_size;
  if (decalage >= 0) {  // this plus grand
    for (int i = substracted.tab_size - 1; i >= 0; i--) {
      temp = this->tab[i + decalage] - substracted.tab[i] - carry;
      carry = this->tab[i + decalage] < (substracted.tab[i] + carry);
      newGrandNombre.tab[i + decalage] = temp % 4294967296;
    }
    // recopie les poids forts de this en propageant la retenue
    for (int i = decalage - 1; i >= 0; i--) {
      temp = this->tab[i] - carry;
      carry = this->tab[i] < carry;
      newGrandNombre.tab[i] = temp % 4294967296;
    }
  } else {  // this plus petit
    for (int i = this->tab_size - 1; i >= 0; i--) {
      temp = this->tab[i] - substracted.tab[i - decalage] - carry;
      carry = this->tab[i] < (substracted.tab[i - decalage] + carry);
      newGrandNombre.tab[i - decalage] = temp % 4294967296;
    }
    // recopie les poids forts de substracted en propageant la retenue
    for (int i = -decalage - 1; i >= 0; i--) {
      temp = substracted.tab[i] - carry;
      carry = substracted.tab[i] < carry;
      newGrandNombre.tab[i] = temp % 4294967296;
    }
  }
  return newGrandNombre;
}

GrandNombre GrandNombre::operator*(const GrandNombre &multiplied) const {
  int newsize = this->tab_size + multiplied.tab_size;
  GrandNombre newGrandNombre = GrandNombre(newsize);
  GrandNombre temp = GrandNombre(newsize);
  for (int i = 0; i < this->tab_size; i++) {
    for (int j = 0; j < multiplied.tab_size; j++) {
      temp.tab[i + j + 1] = this->tab[i] * multiplied.tab[j];
    }
    newGrandNombre = newGrandNombre + temp;
    for (int k = 0; k < temp.tab_size; k++) {
      temp.tab[k] = 0;  // réinitialisation de la variable temp
    }
  }
  return newGrandNombre;
}

bool GrandNombre::operator<=(const GrandNombre &compared) const {
  int i = 0;
  if (this->tab_size >= compared.tab_size) {  // si this plus grand
    int decalage = this->tab_size - compared.tab_size;
    for (int j = 0; j < decalage; j++) {  // là où il n'y a que this
      if (this->tab[j] > 0) {
        return false;
      }
    }
    while (i < compared.tab_size) {
      if (this->tab[i + decalage] < compared.tab[i]) {
        return true;
      } else if (this->tab[i + decalage] > compared.tab[i]) {
        return false;
      }
      i++;
    }
    return true;
  } else {  // si this plus petit
    int decalage = compared.tab_size - this->tab_size;
    for (int j = 0; j < decalage; j++) {  // là où il n'y a que compared
      if (compared.tab[j] > 0) {
        return true;
      }
    }
    while (i < this->tab_size) {
      if (this->tab[i] < compared.tab[i + decalage]) {
        return true;
      } else if (this->tab[i] > compared.tab[i + decalage]) {
        return false;
      }
      i++;
    }
    return true;
  }
}

GrandNombre GrandNombre::operator<<(int decalage) const {
  /* Attention, decalage <= 32 */
  GrandNombre newGrandNombre = GrandNombre(this->tab_size);
  unsigned long long int temp;
  unsigned long long int carry = 0;
  for (int i = this->tab_size - 1; i >= 0; i--) {
    temp = (this->tab[i] << decalage) + carry;
    newGrandNombre.tab[i] = temp % 4294967296;
    carry = temp / 4294967296;
  }
  return newGrandNombre;
}

GrandNombre GrandNombre::operator>>(int decalage) const {
  /* Attention, decalage <= 32 */
  GrandNombre newGrandNombre = GrandNombre(this->tab_size);
  unsigned long long int temp;
  unsigned long long int carry = 0;
  for (int i = 0; i < this->tab_size; i++) {
    temp = (this->tab[i] >> decalage) + (carry << (32 - decalage));
    newGrandNombre.tab[i] = temp % 4294967296;
    carry = this->tab[i] % (1 << decalage);
  }
  return newGrandNombre;
}

GrandNombre GrandNombre::operator%(const GrandNombre &modulo) const {
  GrandNombre newGrandNombre = GrandNombre(modulo.tab_size);
  GrandNombre temp = GrandNombre(modulo.tab_size + 1);
  for (int i = 0; i < this->tab_size; i++) {
    for (int j = 32; j >= 1; j--) {
      temp = temp << 1;
      temp.tab[temp.tab_size - 1] +=
          (this->tab[i] >> (j - 1)) %
          2;  // ajoute le bit 2^j-1 de la case i sur le poids faible de temp
      if (modulo <= temp) {  // retire 2^(tab_size*32 - i) * modulo
        temp = temp - modulo;
      }
    }
  }
  // recopie temp
  for (int i = 0; i < newGrandNombre.tab_size; i++) {
    newGrandNombre.tab[i] = temp.tab[i + 1];
  }
  return newGrandNombre;
}

GrandNombre GrandNombre::subMod(const GrandNombre &substracted,
                                const GrandNombre &modulo) const {
  GrandNombre newGrandNombre = GrandNombre(modulo.tab_size);
  if (substracted <= *this) {
    newGrandNombre = *this - substracted;
  } else {  // réalise soustraction positive
    newGrandNombre = (*this + modulo) - substracted;
  }
  return newGrandNombre;
}

GrandNombre GrandNombre::Montgomery(const GrandNombre &multiplied,
                                    const GrandNombre &R, const GrandNombre &N,
                                    const GrandNombre &V,
                                    const GrandNombre &invV) const {
  GrandNombre T = *this * multiplied;
  GrandNombre S = (T * V) % R;
  GrandNombre M = (S * invV) + T;
  GrandNombre U = GrandNombre(M.tab_size - R.tab_size + 1);
  M = M >> 1;  // la case de R non nulle vaut 2 et pas 1
  for (int i = 0; i < U.tab_size; i++) {
    U.tab[i] = M.tab[i];
  }
  if (N <= U) {
    U = U - N;
  }
  return U;
}

GrandNombre GrandNombre::Phi(const GrandNombre &R, const GrandNombre &N) const {
  return (*this * R) % N;
}

GrandNombre GrandNombre::Pow(const GrandNombre &power, const GrandNombre &R,
                             const GrandNombre &N, const GrandNombre &V,
                             const GrandNombre &invV) const {
  GrandNombre P = R - N;
  for (int i = 0; i < power.tab_size * 32; i++) {
    P = P.Montgomery(P, R, N, V, invV);
    if ((power.tab[i / 32] >> (31 - (i % 32))) % 2 == 1) {
      P = P.Montgomery(*this, R, N, V, invV);
    }
  }
  return P;
}

int main(void) {
  unsigned long long int P[17] = {1, 0, 0, 0, 0, 0, 0, 0,   0,
                                  0, 0, 0, 0, 0, 0, 0, 0x4B};
  unsigned long long int Q[17] = {1, 0, 0, 0, 0, 0, 0, 0,   0,
                                  0, 0, 0, 0, 0, 0, 0, 0x91};
  unsigned long long int E[4] = {12, 2670501072, 1182068202, 1073741881};
  unsigned long long int d_tab[32] = {
      0x0E73DAA8, 0xED437460, 0x3D9F7D9F, 0x3441A270, 0xC58E3328, 0x22920E15,
      0xBB631C02, 0x2B591A14, 0x16DC833D, 0x4F81FFEC, 0xD859A319, 0x9F6DBCB1,
      0xEC84691E, 0x307ED073, 0xFC0CE6CF, 0xF369DDC4, 0xD78AED07, 0x3A37F128,
      0xA0ED18F3, 0x58815C36, 0xBCD6107E, 0x8D0EC25D, 0x2811355D, 0x7DBC0362,
      0x2502FE81, 0x347C5175, 0x2B7C9613, 0x34A05B19, 0xA29C7148, 0x14CE0B0C,
      0x992A46E1, 0xAAB46DC9};
  unsigned long long int r[33] = {2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  unsigned long long int v[32] = {
      0x21D9FC9C, 0x5725A6EC, 0x1D26D8B9, 0xCF136681, 0x15B0C25B, 0xD3D8B3DD,
      0x9F521396, 0x57305FEC, 0x170AD5C3, 0x3E4AAB8C, 0x22D9CAF2, 0x18BEF909,
      0xB6E92D11, 0xD4723A4E, 0x659D4A66, 0xD8FC4E54, 0xBDA3DD10, 0xF8F622FB,
      0x1480E803, 0x33939EEB, 0x8B2573D7, 0x6AA928FB, 0xC94AE1A8, 0x19CCD2D3,
      0x88524BF0, 0xD710442C, 0x29562E47, 0xB41B3666, 0xC6D21EBE, 0x58A70A9D,
      0xDC5C2EFC, 0x7E14DB4D};
  unsigned long long int m[33] = {0, 2147483648, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                  0, 0,          0, 0, 0, 0, 0, 0, 0, 0, 0,
                                  0, 0,          0, 0, 0, 0, 0, 0, 0, 0, 0};
  unsigned long long int un[1] = {1};

  GrandNombre p = GrandNombre(17, P);
  GrandNombre q = GrandNombre(17, Q);
  GrandNombre e = GrandNombre(4, E);
  GrandNombre d = GrandNombre(32, d_tab);
  GrandNombre R = GrandNombre(33, r);
  GrandNombre V = GrandNombre(32, v);
  GrandNombre M = GrandNombre(33, m);
  GrandNombre UN = GrandNombre(1, un);
  GrandNombre N = GrandNombre(33);
  GrandNombre invV = GrandNombre(32);
  GrandNombre N1 = GrandNombre(34);
  GrandNombre invV1 = GrandNombre(33);

  // cout << "variables" << endl;
  // p.Print();
  // q.Print();
  // e.Print();
  // d.Print();
  // R.Print();
  // V.Print();

  N1 = p * q;
  N = N1.resize(33);  // N sur 1025 bits
  invV1 = R - N;
  invV = invV1.resize(32);  // invV < 2^1024
  // N.Print();
  // invV.Print();
  // cout << endl;

  /* exécution */

  // cout << "M : ";
  // M.Print();

  // chrono::steady_clock::time_point begin = chrono::steady_clock::now();
  // M = M.Phi(R, N);
  // GrandNombre C = M.Pow(e, R, N, V, invV);
  // C = C.Montgomery(UN, R, N, V, invV);  // phi-1
  // chrono::steady_clock::time_point end = chrono::steady_clock::now();

  // cout << "C : ";
  // C.Print();

  // chrono::steady_clock::time_point begin2 = chrono::steady_clock::now();
  // C = C.Phi(R, N);
  // GrandNombre D = C.Pow(d, R, N, V, invV);
  // D = D.Montgomery(UN, R, N, V, invV);  // phi-1
  // chrono::steady_clock::time_point end2 = chrono::steady_clock::now();

  // cout << "D : ";
  // D.Print();

  // cout << "Chiffrement = "
  //      << chrono::duration_cast<chrono::milliseconds>(end - begin).count()
  //      << "[ms]" << endl;
  // cout << "Déchiffrement = "
  //      << chrono::duration_cast<chrono::milliseconds>(end2 - begin2).count()
  //      << "[ms]" << endl;

  /* tests */

  // unsigned long long int TEST[16] = {0, 0xFFFFFFFF, 0, 0, 0, 0, 0, 0,
  //                                    0, 0,          0, 0, 0, 0, 2, 1};
  // GrandNombre test = GrandNombre(16, TEST);
  // cout << endl;
  // test.Print();
  // cout << "addition" << endl;
  // GrandNombre n1 = test + test;
  // n1.Print();
  // cout << "soustractions non modulaire" << endl;
  // GrandNombre n2 = p - q;
  // GrandNombre n3 = q - p;
  // n2.Print();
  // n3.Print();
  // GrandNombre n4 = n2 + n3;
  // n4.Print();  // théoriquement 0
  // cout << "multiplication" << endl;
  // GrandNombre n5 = test * test;
  // n5.Print();
  // cout << "decalage" << endl;
  // GrandNombre n8 = test << 2;
  // n8.Print();
  // (test >> 1).Print();
  // (n2 << 1).Print();  // décalage de la soustraction négative
  // cout << "modulo" << endl;
  // GrandNombre n6 = (p + q) % p;
  // GrandNombre n7 = (p + q) % q;
  // n6.Print();
  // n7.Print();
  // cout << "Montgomery" << endl;
  // test.Montgomery(test, R, N, V, invV).Print();

  /* performances */
  cout << endl << "Performances" << endl;

  // Seed
  std::random_device rd;
  // Random number generator
  std::default_random_engine generator(rd());
  // Distribution on which to apply the generator
  std::uniform_int_distribution<long long unsigned> distribution(0, 0xFFFFFFFF);

  // génération des nombres aléatoires
  GrandNombre Liste[100];
  GrandNombre C, D;
  unsigned long long int tab[33];
  for (int i = 0; i < 100; i++) {
    for (int j = 0; j < 32; j++) {
      tab[j + 1] = distribution(generator);
    }
    Liste[i] = GrandNombre(33, tab);
  }

  chrono::steady_clock::time_point begin3 = chrono::steady_clock::now();
  for (int i = 0; i < 100; i++) {
    M = Liste[i].Phi(R, N);
    C = M.Pow(e, R, N, V, invV);
    D = C.Pow(d, R, N, V, invV);
    D = D.Montgomery(UN, R, N, V, invV);
  }
  chrono::steady_clock::time_point end3 = chrono::steady_clock::now();

  cout << "100 RSA = "
       << chrono::duration_cast<chrono::nanoseconds>(end3 - begin3).count()
       << "[ns]" << endl;

  return 0;
}