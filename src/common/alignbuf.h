#ifndef ALIGNBUF_H
#define ALIGNBUF_H

template <typename T, std::size_t align_t=64>
struct align_alloc {
    using value_type = T;

    template <typename U>
    struct rebind {
        using other = align_alloc<U, align_t>;
    };

    constexpr align_alloc() noexcept = default;
    constexpr align_alloc(const align_alloc &) noexcept = default;


    template <typename U>
    constexpr align_alloc(const align_alloc<U, align_t> &) noexcept {}

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

#endif
