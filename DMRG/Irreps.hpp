#pragma once

#include <string>

class Irreps {
   private:
    bool isActivated;
    int groupNumber;
    int irrepsNumber;

    void init(const int group);

    static std::string getGroupName(const int group);
    static std::string getIrrepName(const int group, const int irrep);

   public:
    Irreps();
    Irreps(const int group);
    virtual ~Irreps();
    
    bool setGroup(const int group);
    bool getIsActivated() const;
    
    int getGroupNumber() const;
    int numberOfIrreps() const;
    
    std::string groupName() const;
    std::string irrepName(const int irrep) const;
    
    void symmFillArray(int* const arr) const;
    
    static int numberOfIrreps(const int group);
    static int directProduct(const int irrepFirst, const int irrepSecond);
    static std::string groupName(const int group);
    static void symmFillArray(int* const arr, const std::string symmLabel);
    static void printAll();
};
