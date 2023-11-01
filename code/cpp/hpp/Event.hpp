

template <typename ing, typename num>
struct Event
{
  ing id; // Event index, not necessarily the true event ID.
  ing size;
  ing *ind;
  num *val;
  double sumOfSquaredVal;
  Event (){ sumOfSquaredVal = -1; }
  void reset(ing *ind, num *val, ing id, ing size)
  {
    this->id = id;
    this->ind = ind;
    this->val = val;
    this->size = size;
    sumOfSquaredVal = -1;
  }
  Event( ing *ind, num *val, ing id, ing size) { reset(ind, val, id, size); }
};














