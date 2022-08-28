#ifndef SOCKET_WRAPPER
#define SOCKET_WRAPPER

#include <netdb.h>
#include <netinet/in.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <unistd.h>

#include <cstring>
#include <exception>
#include <string>

class SocketException : public std::exception {
    const char* what() const throw() {
        return "Socket error";
    }
};

class Socket {
   public:
    Socket(int fd);
    Socket();
    ~Socket();
    void send(std::string msg);
    std::string recv(uint16_t len);

   protected:
    int socketFd;
};

class ClientSocket : public Socket {
   public:
    ClientSocket(std::string hostname, uint16_t port);
};

class ServerSocket : public Socket {
   public:
    ServerSocket(std::string hostname, uint16_t port);
    Socket* acceptClient();
};

#endif