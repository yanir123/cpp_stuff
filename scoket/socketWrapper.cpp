#include "socketWrapper.hpp"

Socket::Socket(int fd) {
    this->socketFd = fd;
}

Socket::~Socket() {
    if (this->socketFd >= 0) {
        close(this->socketFd);
    }
}

ClientSocket::ClientSocket(std::string hostname, uint16_t port) {
    struct hostent* server;
    struct sockaddr_in serv_addr;

    this->socketFd = socket(AF_INET, SOCK_STREAM, 0);

    if (this->socketFd < 0) {
        throw SocketException();
    }

    server = gethostbyname(hostname.c_str());
    if (!server) {
        throw SocketException();
    }

    bzero(&serv_addr, sizeof(serv_addr));
    serv_addr.sin_family = AF_INET;
    bcopy(server->h_addr, &serv_addr.sin_addr.s_addr, server->h_length);
    serv_addr.sin_port = htons(port);

    if (connect(this->socketFd, (struct sockaddr*)&serv_addr, sizeof(serv_addr)) < 0) {
        throw SocketException();
    }
}

void Socket::send(std::string msg) {
    ssize_t bytes_written = 0;
    bytes_written = write(this->socketFd, msg.c_str(), msg.length());
    if (bytes_written < 0) {
        throw SocketException();
    }
}

std::string Socket::recv(uint16_t len) {
    std::string res = "";
    ssize_t bytes_read = 0;
    char* buffer = new char[len]();

    bytes_read = read(this->socketFd, buffer, len);
    if (bytes_read < 0) {
        throw SocketException();
    }
    res += buffer;

    delete buffer;

    return res;
}

ClientSocket::~ClientSocket() {
    close(this->socketFd);
}

ServerSocket::ServerSocket(std::string hostname, uint16_t port) {
    struct sockaddr_in serv_addr;

    this->socketFd = socket(AF_INET, SOCK_STREAM, 0);
    if (this->socketFd < 0) {
        throw SocketException();
    }

    bzero((char*)&serv_addr, sizeof(serv_addr));
    serv_addr.sin_family = AF_INET;
    serv_addr.sin_addr.s_addr = INADDR_ANY;
    serv_addr.sin_port = htons(port);

    if (bind(this->socketFd, (struct sockaddr*)&serv_addr, sizeof(serv_addr)) < 0) {
        throw SocketException();
    }

    listen(this->socketFd, 5);
}

Socket ServerSocket::acceptClient() {
    struct sockaddr_in cli_addr;
    socklen_t clilen;

    clilen = sizeof(cli_addr);
    int clientScokFd = accept(this->socketFd, (struct sockaddr*)&cli_addr, &clilen);
    if (clientScokFd < 0) {
        throw SocketException();
    }

    return Socket(clientScokFd);
}