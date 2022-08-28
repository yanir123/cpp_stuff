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
/*
    Exception class to alert of network failure
*/
class SocketException : public std::exception {
    const char* what() const throw() {
        return "Socket error";
    }
};

/*
    Class to wrap the basic linux socket interface
*/
class Socket {
   public:
    /*
    Constructor for the socket class with existing socker fd
        @param fd int: open socket file descriptor
    */
    Socket(int fd);

    /*
    Constructor for the socket class with empty initialization
    */
    Socket();

    /*
    Destructor for the socket class
    */
    ~Socket();

    /*
    Socket class function to send data
        @param msg string - message to send over the socket
    */
    void send(std::string msg);

    /*
    Socket function to recieve data 
        @param len uint16_t - max length of the expected message
        @return string - message recieved over the socket
    */
    std::string recv(uint16_t len);

   protected:
    int socketFd;
};

/*
Class to represent client side socket
*/
class ClientSocket : public Socket {
   public:
    /*
    Constructor to initialize the socket
        @param hostname string - hostname to connect to
        @param port uint16_t - port to connect on
    */
    ClientSocket(std::string hostname, uint16_t port);
};

/*
Class to represent Server side socket
*/
class ServerSocket : public Socket {
   public:
    /*
    Constructor to initialize the socket
        @param hostname string - hostname to bind to
        @param port uint16_t - port to bind on
    */
    ServerSocket(std::string hostname, uint16_t port);

    /*
    Function to accept clients
        @return Socket - new client
    */
    Socket* acceptClient();
};

#endif