// A 2D/3D collision grid. zlib/libpng licensed
// - rlyeh, 2012-2013

#pragma once

#include <cassert>
#include <map>
#include <set>
#include <tuple>

template< typename T >
class grid3
{
    typedef std::tuple<int,int,int> key;
    typedef typename std::set<T>    set;

    const int cell_size;
    std::map< key, set > map;
    typename std::map< key, set >::iterator cursor;

    inline key at( int x, int y, int z ) const {
        return key( x / cell_size, y / cell_size, z / cell_size );
    }

    inline set &with( int x, int y, int z ) {
        auto idx = at(x,y,z);
        return ( map[ idx ] = map[ idx ] );
    }

    public:

    grid3( int _cell_size = 1 ) : cursor(map.end()), cell_size(_cell_size) {
        assert( _cell_size > 0 );
    }

    void setup( int _cell_size ) {
        *this = grid3( _cell_size );
    }

    bool find( int x, int y = 0, int z = 0 ) {
        cursor = map.find( at(x,y,z) );
        return( cursor != map.end() );
    }

    set &found() const {
        return cursor->second;
    }

    size_t size( int x, int y = 0, int z = 0 ) const {
        return find(x,y,z) ? found().size() : 0;
    }

    size_t count( const T &elem, int x, int y = 0, int z = 0 ) const {
        return find(x,y,z) ? found().count( elem ) : 0;
    }

    void insert( const T &elem, int x, int y = 0, int z = 0 ) {
        with(x,y,z).insert( elem );
    }

    void erase( const T &elem, int x, int y = 0, int z = 0 ) {
        with(x,y,z).erase( elem );
    }

    void clear( int x, int y = 0, int z = 0 ) {
        with(x,y,z) = set(); //.clear();
    }
};
