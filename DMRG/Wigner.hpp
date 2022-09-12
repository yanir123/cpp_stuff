#pragma once

constexpr int WIGNER_FACTORIAL_MAX = 191;
constexpr int WIGNER_MAX_2J = 95;

class Wigner {
   private:
    static const double sqrtFact[WIGNER_FACTORIAL_MAX + 1];
    static bool triangleFails(const int ja, const int jb, const int jc);
    static bool sqrtDelta(const int ja, const int jb, const int jc);

   public:
    Wigner();
    ~Wigner();
    static int maxTwoJ();
    static double wignerThreeJ(const int ja, const int jb, const int jc, const int ma, const int mb, const int mc);
    static double wignerSixJ(const int ja, const int jb, const int jc, const int jd, const int je, const int jf);
    static double wignerNineJ(const int ja, const int jb, const int jc, const int jd, const int je, const int jf, const int jg, const int jh, const int ji);
};
