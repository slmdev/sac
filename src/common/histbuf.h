#ifndef HISTBUF_H
#define HISTBUF_H

#include "../global.h"

template <typename T, std::size_t align_t>
struct align_alloc {
    using value_type = T;

    template <typename U>
    struct rebind {
        using other = align_alloc<U, align_t>;
    };

    align_alloc() noexcept = default;

    template <typename U>
    constexpr align_alloc(const align_alloc<U, align_t>&) noexcept {}

    T* allocate(std::size_t n) {
        auto ptr = static_cast<T*>(::operator new(n * sizeof(T), std::align_val_t(align_t)));
        return ptr;
    }

    void deallocate(T* p, std::size_t n) noexcept {
        ::operator delete(p, n * sizeof(T), std::align_val_t(align_t));
    }
};

template <typename T, typename U, std::size_t align_t>
bool operator==(const align_alloc<T, align_t>&, const align_alloc<U, align_t>&) noexcept {
    return true;
}

template <typename T, typename U, std::size_t align_t>
bool operator!=(const align_alloc<T, align_t>&, const align_alloc<U, align_t>&) noexcept {
    return false;
}


// rolling buffer, n must be power of two
template <typename T>
class HistBuffer {
  public:
      HistBuffer(int n)
      :n(n),mask(n-1),pos(0),buf(n)
      {
      }
      T& operator[](int idx) { return buf[(pos-idx)&mask];};
      const T operator[](int idx) const { return buf[(pos-idx)&mask];};
      void PushBack(T val)
      {
         buf[pos]=val;
         pos=(pos+1)&mask;
      }
  private:
    int n,mask,pos;
    std::vector <T> buf;
};

// circulating buffer
// operator [] starts from 0=newest to oldest
template <typename T>
class RollBuffer {
  public:
    RollBuffer(int n)
    :n(n),pos(-1),buf(n)
    {

    }
    const T& operator[](int idx) const
    {
      return buf[clamp_idx(pos-idx)];
    };
    void push(T val)
    {
      if (++pos>=n) pos=0;
      buf[pos]=val;
    }
    int clamp_idx(int idx) const
    {
      if (idx>=n) idx-=n;
      else if (idx<0) idx+=n;
      return idx;
    }
    const std::vector<T> &getbuf() {return buf;};
  private:
    int n,pos;
    std::vector <T> buf;
};

//circulating buffer, bi-partit
//operator [] starts from 0=newest to oldest
template <typename T>
class RollBuffer2 {
  public:
    RollBuffer2(std::size_t capacity)
    :n(capacity),pos(0),buf(2*capacity)
    {
      #if 0
        auto* data_ptr = buf.data();
        std::uintptr_t address = reinterpret_cast<std::uintptr_t>(data_ptr);
        std::size_t align = 1ULL << __builtin_ctzll(address);
        std::cout << " rb2 adr: " << static_cast<void*>(data_ptr) << ", align: " << align << '\n';
      #endif
    }
    void push(T val)
    {
      pos = (pos + n - 1) % n;
      buf[pos] = val;
      buf[pos+n] = val;
    }

    const T& operator[](int index) const
    {
      return buf[pos + index];
    }

    const std::span<T> get_span() const {
      return std::span<const>{buf.data() + pos,n};
    }
    const T* data() const {
      return buf.data() + pos;
    }
  private:
    std::size_t n,pos;
    std::vector<T> buf;
};


#endif // HISTBUF_H
