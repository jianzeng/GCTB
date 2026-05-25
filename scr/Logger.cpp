/*
   Global Singleton logger system with robust functions and thread safe.

      * open a logger file first, otherwise no output into log
      * e(int level, string message):  prompt error, and exit the program, level is the number of indent spaces
      * i:  prompt information
      * w:  prompt warning message
      * d:  debug message, only seen in the debug mode
      * m:  message that only show on the terminal that not log into logger file
      * l:  log into logger file only;
      * p:  progress, that always show in one line in the terminal but no output into logger file.
      * << :  use like std::cout

   Developed by Zhili Zheng<zhilizheng@outlook.com>

   This file is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   A copy of the GNU General Public License is attached along with this program.
   If not, see <http://www.gnu.org/licenses/>.
*/

#include "Logger.hpp"
#include <iostream>
#include <cstdlib>
#include <streambuf>

#ifdef _WIN32
#include <io.h>
#include <stdio.h>
#define ISATTY(x) _isatty(_fileno(x))
#else 
#include <unistd.h>
#define ISATTY(x) isatty(fileno(x))
#endif

using std::cout;
using std::cerr;

std::mutex Logger::log_mutex;
std::string Logger::empty = {};
std::map<string, std::chrono::time_point<std::chrono::steady_clock>> Logger::time_map;
#ifdef _WIN32 
std::map<Logger::Type, string> Logger::style_map = {{Logger::INFO, ""}, {Logger::PROMPT, ""},
                                                     {Logger::PROGRESS, "\r"}, {Logger::WARN, ""},
                                                     {Logger::ERROR, ""}, {Logger::DEBUG, ""}};
#else
std::map<Logger::Type, string> Logger::style_map = {{Logger::INFO, "\033[0m"}, {Logger::PROMPT, "\033[0;32m"},
                                                     {Logger::PROGRESS, "\r"}, {Logger::WARN, "\033[0;33m"},
                                                     {Logger::ERROR, "\033[0;31m"}, {Logger::DEBUG, "\033[0;34m"}};
#endif
Logger* Logger::m_pThis = NULL;
Logger::Type Logger::m_stat = Logger::INFO;
std::ofstream Logger::m_logFile;
string Logger::m_FileName = "";
std::streambuf* Logger::s_savedCoutBuf = NULL;
std::streambuf* Logger::s_savedCerrBuf = NULL;
std::streambuf* Logger::s_coutTeeBuf = NULL;
std::streambuf* Logger::s_cerrTeeBuf = NULL;

namespace {

/// Copy every character from cout/cerr to the log file (including empty lines).
class StreamTeeBuf : public std::streambuf {
    std::streambuf* dest_;
    std::ofstream* log_;
public:
    StreamTeeBuf(std::streambuf* dest, std::ofstream* log) : dest_(dest), log_(log) {}
protected:
    int overflow(int c) {
        if (c == traits_type::eof())
            return traits_type::not_eof(c);
        const char ch = traits_type::to_char_type(c);
        if (dest_) dest_->sputc(ch);
        if (log_ && log_->is_open()) log_->put(ch);
        return c;
    }
    int sync() {
        if (dest_) dest_->pubsync();
        if (log_ && log_->is_open()) log_->flush();
        return 0;
    }
};

} // namespace

Logger::Logger(){
    if(!ISATTY(stdout)){
        Logger::style_map =  {{Logger::INFO, ""}, {Logger::PROMPT, ""},
                              {Logger::PROGRESS, "\n"}, {Logger::WARN, ""},
                              {Logger::ERROR, ""}, {Logger::DEBUG, ""}};
    }

}

Logger* Logger::GetLogger(){
    if(m_pThis == NULL) {
        m_pThis = new Logger();
    }
    return m_pThis;
}

bool Logger::check(){
    if(m_pThis == NULL){
        return false;
    }
    if(!m_logFile.is_open()){
        return false;
    }
    if(m_FileName.empty()){
        return false;
    }
    return true;
}

void Logger::open(string ofile){
    if(check()){
        cout << "Logger has been set, not support to set another time" << endl;
    }else{
        m_logFile.open(ofile, std::ios::out);
        m_FileName = ofile;
        if(!check()){
            m_pThis->e(0, "can't write to log file [" + ofile + "].\nPlease check file/folder permission or disk quota.");
        }
    }
}

void Logger::attachStdout(){
    if (!check() || s_savedCoutBuf != NULL)
        return;
    GetLogger();
    s_savedCoutBuf = std::cout.rdbuf();
    s_coutTeeBuf = new StreamTeeBuf(s_savedCoutBuf, &m_logFile);
    std::cout.rdbuf(s_coutTeeBuf);
    s_savedCerrBuf = std::cerr.rdbuf();
    s_cerrTeeBuf = new StreamTeeBuf(s_savedCerrBuf, &m_logFile);
    std::cerr.rdbuf(s_cerrTeeBuf);
}

bool Logger::stdoutIsAttached(){
    return s_savedCoutBuf != NULL;
}

void Logger::restoreStdout(){
    if (s_coutTeeBuf != NULL)
        std::cout.flush();
    if (s_cerrTeeBuf != NULL)
        std::cerr.flush();
    if (s_savedCoutBuf != NULL) {
        std::cout.rdbuf(s_savedCoutBuf);
        s_savedCoutBuf = NULL;
    }
    delete s_coutTeeBuf;
    s_coutTeeBuf = NULL;
    if (s_savedCerrBuf != NULL) {
        std::cerr.rdbuf(s_savedCerrBuf);
        s_savedCerrBuf = NULL;
    }
    delete s_cerrTeeBuf;
    s_cerrTeeBuf = NULL;
}

void Logger::close(){
    restoreStdout();
    if (m_logFile.is_open())
        m_logFile.close();
    m_FileName.clear();
}

void Logger::flush(){
    m_logFile.flush();
}

void Logger::LogUnlocked(int level, Type type, const string& prompt, const string& message){
    string spaces(level * 2, ' ');
    (*m_pThis) << spaces << type << prompt << INFO << message << endl;
}

void Logger::Log(int level, Type type, const string& prompt, const string& message){
    std::lock_guard<std::mutex> lock(log_mutex);
    LogUnlocked(level, type, prompt, message);
}

void Logger::e(int level, const string& message, const string& title){
    string head = title.empty() ? "Error: " : (title+" ");
    std::lock_guard<std::mutex> lock(log_mutex);
    if (stdoutIsAttached()) {
        cout << " \n" << message << endl;
        cout.flush();
        m_logFile.flush();
    } else {
        LogUnlocked(level, ERROR, head, message);
        LogUnlocked(level, INFO, "", "An error occurs, please check the options or data");
    }
    exit(EXIT_FAILURE);
}

int Logger::precision(int p){
    cout.precision(p);
    if (!stdoutIsAttached() && m_logFile.is_open())
        m_logFile.precision(p);
    return cout.precision();
}

string Logger::setprecision(int p){
    cout.precision(p);
    if (!stdoutIsAttached() && m_logFile.is_open())
        m_logFile.precision(p);
    return("");
}

int Logger::precision(){
    return cout.precision();
}

void Logger::i(int level, const string& message, const string& title){
    string head = title.empty() ? "" : (title + " ");
    std::lock_guard<std::mutex> lock(log_mutex);
    if (stdoutIsAttached()) {
        cout << head << message << endl;
    } else {
        LogUnlocked(level, PROMPT, head, message);
    }
}

void Logger::w(int level, const string &message, const string& title) {
    string head = title.empty() ? "Warning: " : (title + " ");
    std::lock_guard<std::mutex> lock(log_mutex);
    if (stdoutIsAttached()) {
        cerr << head << message << endl;
    } else {
        LogUnlocked(level, WARN, head, message);
    }
}

void Logger::d(int level, const string &message, const string& title) {
    #ifndef NDEBUG
    string head = title.empty() ? "Debug: " : (title + " ");
    m_pThis->Log(level, DEBUG, head, message);
    #endif
}

void Logger::p(int level, const string &message, const string& title) {
    string head = title.empty() ? "" : (title + " ");
    string spaces(level * 2, ' ');
    m_stat = PROGRESS;
    std::lock_guard<std::mutex> lock(log_mutex);
    (*m_pThis) << spaces << head << message << PROGRESS << std::flush;
    m_stat = INFO;
}

void Logger::m(int level, const string &message, const string& title){
    string head = title.empty() ? "" : (title + " ");
    string spaces(level * 2, ' ');
    m_stat = PROGRESS;
    std::lock_guard<std::mutex> lock(log_mutex);
    (*m_pThis) << spaces << PROMPT << head << INFO << message << endl;
    m_stat = INFO;
}

void Logger::l(int level, const string &message, const string &title){
    string head = title.empty() ? "" : (title + " ");
    string spaces(level * 2, ' ');
    m_logFile << spaces << head << message << endl;
}

void Logger::ts(string key){
    time_map[key] = std::chrono::steady_clock::now();
}

float Logger::tp(string key){
    auto end = std::chrono::steady_clock::now();
    float secs = 0.0;
    if(time_map.find(key) != time_map.end()){
        auto duration = end - time_map[key];
        secs = std::chrono::duration_cast<std::chrono::duration<float>>(duration).count();
    }else{
        m_pThis->w(2, "can't find key of time start");
    }
    return secs;
}

Logger& Logger::operator<<(Type type){
    m_stat = type;
    cout << style_map[type];
    return *m_pThis;
}

Logger& Logger::operator<<(std::ostream& (*op)(std::ostream&)){
    (*op)(cout);
    if(m_stat != PROGRESS && m_logFile.is_open() && !stdoutIsAttached()){
        (*op)(m_logFile);
    }
    return *m_pThis;
}

Logger& Logger::operator<<(std::ios& (*pf)(std::ios&)){
    cout << pf;
    if (m_logFile.is_open() && !stdoutIsAttached())
        m_logFile << pf;
    return *m_pThis;
}

Logger& Logger::operator<<(std::ios_base& (*pf)(std::ios_base&)){
    cout << pf;
    if (m_logFile.is_open() && !stdoutIsAttached())
        m_logFile << pf;
    return *m_pThis;
}


