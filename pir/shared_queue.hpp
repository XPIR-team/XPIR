/* Copyright (C) 2014 Carlos Aguilar Melchor, Joris Barrier, Marc-Olivier Killijian
 * This file is part of XPIR.
 *
 *  XPIR is free software: you can redistribute it and/or modify
 *	it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  XPIR is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with XPIR.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef DEF_SHARED_BUFFER
#define DEF_SHARED_BUFFER

#include <stdio.h>
#include <stdlib.h>
#include <limits>
#include <string>
#include <boost/thread.hpp>

#ifdef __APPLE__
#include <boost/interprocess/sync/named_semaphore.hpp>
#else
#include <boost/interprocess/sync/interprocess_semaphore.hpp>
#endif
#include <boost/pending/queue.hpp>
using namespace boost;
using namespace boost::interprocess;
using namespace std;

template<typename T> 
class shared_queue
{

  public:
    shared_queue(const std::string& name, unsigned int max_size=SEM_VALUE_MAX);
    ~shared_queue();
    std::string id;

    void push(T);
    void pop();
    unsigned int size();
    T front();
    T pop_front();
    bool empty();
  private:
		boost::mutex mutex;
#ifdef __APPLE__
    named_semaphore num_stored, num_space;
#else
    interprocess_semaphore num_stored, num_space;
#endif
    //Items to fill
    boost::queue<T> data;
};

template <typename T>
shared_queue<T>::shared_queue(const std::string& name, unsigned int max_size):
#ifdef __APPLE__
  num_stored(open_or_create, string(name + "_num_stored").c_str(), 0, permissions(777)), 
  num_space(open_or_create,  string(name + "_num_space").c_str(), max_size, permissions(777)),
  id(name)
#else
 num_stored(0), num_space(max_size)
#endif
{
#ifdef __APPLE__
  /*reset semaphore*/
  while (num_stored.try_wait()){}

  if (num_space.try_wait() == false) {
    for (unsigned int i = 0 ; i < max_size ; i++) {
      num_space.post();
    }
  }
  else {
      num_space.post();
    }
#endif
}

// template <typename T>
// shared_queue<T>::shared_queue(const std::string& name): 
//   shared_queue(name, SEM_VALUE_MAX)
//   {}

template <typename T>
void shared_queue<T>::push(T item){
  num_space.wait();
  mutex.lock();
  data.push(item);
  mutex.unlock();
  num_stored.post();
}

template <typename T>
T shared_queue<T>::front(){
  T pt;
  num_stored.wait();
  mutex.lock();
  pt = data.front();
  mutex.unlock();
  num_stored.post();
  return pt;
}

template <typename T>
T shared_queue<T>::pop_front(){
  T pt;
  num_stored.wait();
  mutex.lock();
  pt = data.front();
  data.pop();
  mutex.unlock();
  num_space.post();
  return pt;
}

template <typename T>
void shared_queue<T>::pop(){
  num_stored.wait();
  mutex.lock();
  data.pop();
  mutex.unlock();
  num_space.post();
}

template <typename T>
bool shared_queue<T>::empty(){
  bool b;
  mutex.lock();
  b = data.empty();
  mutex.unlock();
  return b;
}

template <typename T>
unsigned int shared_queue<T>::size(){
  int s;
  mutex.lock();
  s = data.size();
  mutex.unlock();
  return s;
}

template <typename T>
shared_queue<T>::~shared_queue(){
#ifdef __APPLE__
  boost::interprocess::named_semaphore::remove(string(id + "_mutex").c_str());
  boost::interprocess::named_semaphore::remove(string(id + "_num_stored").c_str());
  boost::interprocess::named_semaphore::remove(string(id + "_num_space").c_str());
#endif
}

#endif
