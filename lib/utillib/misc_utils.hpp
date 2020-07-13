#ifndef HAVE_MISC_UTILS_HPP
#define HAVE_MISC_UTILS_HPP
#endif



template <class XX> int delta(XX firstarg,XX secondarg){

  /* implementation of Kroenecker-Delta for any type object
   * == operator must be properly overloader for that object type
   * Nachiket Gokhale gokhalen@bu.edu Mar 8 2005
   * This should be made a more elegant oneliner.
   */

  if(firstarg==secondarg){
    return 1;
  }
  else{
    return 0;
  }

  return 0; // never get here, this line exists to trick compiler  into thinking that the function closes with a return

}
