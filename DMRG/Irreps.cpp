#include "Irreps.hpp"

#include <cmath>
#include <iostream>

Irreps::Irreps() {
    this->isActivated = false;
}

Irreps::Irreps(const int group) {
    this->init(group);
}

Irreps::~Irreps() {
}

void Irreps::init(const int group) {
    if (group >= 0 && group <= 7) {
        this->isActivated = true;
        this->groupNumber = group;
        this->irrepsNumber = this->numberOfIrreps(group);
    } else {
        this->isActivated = false;
    }
}

bool Irreps::setGroup(const int group) {
    this->init(group);

    return this->isActivated;
}

bool Irreps::getIsActivated() const {
    return this->isActivated;
}

int Irreps::getGroupNumber() const {
    return this->isActivated ? this->groupNumber : -1;
}

std::string Irreps::groupName() {
    return this->isActivated ? getGroupName(this->groupNumber) : "error";
}

std::string Irreps::groupName(const int group) {
    return (group >= 0 && group <= 7) ? getGroupName(group) : "error";
}

std::string Irreps::getGroupName(const int group) {
    std::string names[8] = {"c1", "ci", "c2", "cs", "d2", "c2v", "c2h", "d2h"};

    if (group < 0 || group > 7) {
        return "error";
    }

    return names[group];
}

int Irreps::numberOfIrreps() const {
    return this->isActivated ? this->irrepsNumber : -1;
}

int Irreps::numberOfIrreps(const int group) {
    if (group < 0 || group > 7) {
        return -1;
    }

    return (int)std::pow(2, std::ceil((float)group / 3));
}

std::string Irreps::irrepName(const int irrep) const {
    if (!this->isActivated) {
        return "error1";
    }

    if (irrep < 0 || irrep >= this->irrepsNumber) {
        return "error2";
    }

    return getIrrepName(this->groupNumber, irrep);
}

std::string Irreps::getIrrepName(const int group, const int irrep) {
    if (group == 0) {
        if (irrep == 0) {
            return "A";
        }
    }

    if (group == 1) {
        if (irrep == 0) {
            return "Ag";
        }
        if (irrep == 1) {
            return "Au";
        }
    }

    if (group == 2) {
        if (irrep == 0) {
            return "A";
        }
        if (irrep == 1) {
            return "B";
        }
    }

    if (group == 3) {
        if (irrep == 0) {
            return "Ap";
        }
        if (irrep == 1) {
            return "App";
        }
    }

    if (group == 4) {
        if (irrep == 0) {
            return "A";
        }
        if (irrep == 1) {
            return "B1";
        }
        if (irrep == 2) {
            return "B2";
        }
        if (irrep == 3) {
            return "B3";
        }
    }

    if (group == 5) {
        if (irrep == 0) {
            return "A1";
        }
        if (irrep == 1) {
            return "A2";
        }
        if (irrep == 2) {
            return "B1";
        }
        if (irrep == 3) {
            return "B2";
        }
    }

    if (group == 6) {
        if (irrep == 0) {
            return "Ag";
        }
        if (irrep == 1) {
            return "Bg";
        }
        if (irrep == 2) {
            return "Au";
        }
        if (irrep == 3) {
            return "Bu";
        }
    }

    if (group == 7) {
        if (irrep == 0) {
            return "Ag";
        }
        if (irrep == 1) {
            return "B1g";
        }
        if (irrep == 2) {
            return "B2g";
        }
        if (irrep == 3) {
            return "B3g";
        }
        if (irrep == 4) {
            return "Au";
        }
        if (irrep == 5) {
            return "B1u";
        }
        if (irrep == 6) {
            return "B2u";
        }
        if (irrep == 7) {
            return "B3u";
        }
    }

    return "error2";
}

void Irreps::symmFillArray(int* const arr) const {
    if (this->isActivated) {
        symmFillArray(arr, this->groupName());
    }
}

void Irreps::symmFillArray(int* const arr, const std::string symmLabel) {
    if (symmLabel.compare("c1") == 0) {
        arr[0] = 1;
    }
    if ((symmLabel.compare("ci") == 0) || (symmLabel.compare("c2") == 0) || (symmLabel.compare("cs") == 0)) {
        arr[0] = 1;
        arr[1] = 2;
    }
    if ((symmLabel.compare("d2") == 0)) {
        arr[0] = 1;
        arr[1] = 4;
        arr[2] = 3;
        arr[3] = 2;
    }
    if ((symmLabel.compare("c2v") == 0) || (symmLabel.compare("c2h") == 0)) {
        arr[0] = 1;
        arr[1] = 4;
        arr[2] = 2;
        arr[3] = 3;
    }
    if ((symmLabel.compare("d2h") == 0)) {
        arr[0] = 1;
        arr[1] = 4;
        arr[2] = 6;
        arr[3] = 7;
        arr[4] = 8;
        arr[5] = 5;
        arr[6] = 3;
        arr[7] = 2;
    }
}

void Irreps::printAll() {
    for (int thegroup = 0; thegroup < 8; thegroup++) {
        std::cout << "######################################################" << std::endl;
        std::cout << "Name = " << getGroupName(thegroup) << std::endl;
        std::cout << "nIrreps = " << numberOfIrreps(thegroup) << std::endl;
        std::cout << "Multiplication table :" << std::endl;
        for (int irrep1 = -1; irrep1 < numberOfIrreps(thegroup); irrep1++) {
            for (int irrep2 = -1; irrep2 < numberOfIrreps(thegroup); irrep2++) {
                if ((irrep1 == -1) && (irrep2 == -1)) std::cout << "\t";
                if ((irrep1 == -1) && (irrep2 >= 0)) std::cout << getIrrepName(thegroup, irrep2) << "\t";
                if ((irrep2 == -1) && (irrep1 >= 0)) std::cout << getIrrepName(thegroup, irrep1) << "\t";
                if ((irrep2 >= 0) && (irrep1 >= 0)) std::cout << getIrrepName(thegroup, directProduct(irrep1, irrep2)) << "\t";
            }
            std::cout << std::endl;
        }
    }
    std::cout << "######################################################" << std::endl;
}