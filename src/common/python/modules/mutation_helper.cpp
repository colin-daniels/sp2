// FIXME ugh, the intent was to avoid this import, but I didn't manage
//       to break up fake_modules.hpp quite enough for that.
#include "Python.h"

#include "common/python/modules/fake_modules.hpp"

namespace {

const char* name = "mutation_helper";
const char* text = "\n\
from collections import defaultdict                                                        \n\
                                                                                           \n\
class supercell_index_mapper:                                                              \n\
    # Takes:                                                                               \n\
    #   reverse_image_map:  list with an integer for each atom in the supercell,           \n\
    #                       indicating which atom in the primary cell it is an image of.   \n\
    def __init__(self, reverse_image_map):                                                 \n\
        images = defaultdict(list)                                                         \n\
        for (image, primary) in enumerate(reverse_image_map):                              \n\
            images[primary].append(image)                                                  \n\
                                                                                           \n\
        self._prim_size = len(images)                                                      \n\
                                                                                           \n\
        self._images = {}                                                                  \n\
        for star in images.values():                                                       \n\
            for x in star:                                                                 \n\
                self._images[x] = star                                                     \n\
                                                                                           \n\
        self._image_rev = reverse_image_map                                                \n\
                                                                                           \n\
        self._unique_prim = sorted(set(self._image_rev))                                   \n\
                                                                                           \n\
    def primitive_size(self):                                                              \n\
        ''' Get the number of primitive cell atoms. '''                                    \n\
        return self._prim_size                                                             \n\
                                                                                           \n\
    def supercell_size(self):                                                              \n\
        ''' Get the number of atoms in the supercell. '''                                  \n\
        return len(self._image_rev)                                                        \n\
                                                                                           \n\
    def images_of(self, idx):                                                              \n\
        '''                                                                                \n\
        Get all indices in the supercell equivalant to an index under                      \n\
        the primitive lattice.                                                             \n\
                                                                                           \n\
        :param idx: primitive cell atom index, or iterable of indices                      \n\
        :return: list of indices                                                           \n\
        '''                                                                                \n\
        if isinstance(idx, int):                                                           \n\
            return self.images_of([idx])                                                   \n\
        elif hasattr(idx, '__iter__'):                                                     \n\
            return [img for i in idx for img in self._images[i]]                           \n\
        else:                                                                              \n\
            raise TypeError('Expected int or iterable of ints')                            \n\
                                                                                           \n\
    def primitive_index(self, idx):                                                        \n\
        ''' Look up primitive cell index for a given atom in the supercell. '''            \n\
        return self._image_rev[idx]                                                        \n\
";

} // anonymous namespace

sp2::python::fake_module_t sp2::python::fake_modules::mutation_helper = {name, text};
