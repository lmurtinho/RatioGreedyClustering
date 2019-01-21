#include <stdio.h>
#include <stdbool.h>


void heap_check_position(HEAP *heap, int i, int n) {
  int p;

  if (i > 0) {
    p = (i - 2 + (i % 2)) / 2;
    if (heap[i].val < heap[p].val) {
      printf("n = %d: Element %d has value %.10e smaller than parent %d %.10e\n",
        n, i, heap[i].val, p, heap[p].val);
    }
  }

  p = 2*i+1;
  if ( (p < n) && (heap[i].val > heap[p].val) ) {
    printf("n =n %d: Element %d has value %.10e larger than son %d %.10e\n",
      n, i, heap[i].val, p, heap[p].val);
  }

  p++;
  if ( (p < n) && (heap[i].val > heap[p].val) ) {
    printf("n = %d: Element %d has value larger than son %d\n", n, i, p);
  }

  if ( (heap[i].prev >= 0) && (heap[heap[i].prev].next != i) ) {
    printf("n = %d: Previous element %d has next %d instead of %d\n",
      n, heap[i].prev, heap[heap[i].prev].next, i);
  }

  if ( (heap[i].next >= 0) && (heap[heap[i].next].prev != i) ) {
    printf("n = %d: Next element %d has prev %d instead of %d\n",
      n, heap[i].next, heap[heap[i].next].prev, i);
  }

  if ( (heap[i].next > (n-1) ) ) {
    printf("Element %d next %d is outside heap with %d elements\n",
      i, heap[i].next, n);
  }

  if ( (heap[i].prev > (n-1) ) ) {
    printf("Element %d prev %d is outside heap with %d elements\n",
      i, heap[i].prev, n);
  }

  if ( (n>1) && (heap[i].prev == heap[i].next) ) {
    printf("n = %d: element %d has same prev %d and next %d\n",
      n, i, heap[i].prev, heap[i].next);
  }

  if (heap[i].ini >= heap[i].end) {
    printf("n = %d: element %d has ini %d >= end %d\n",
      n, i, heap[i].ini, heap[i].end);
  }

  if (heap[i].ini >= heap[i].end) {
    printf("n = %d: element %d has ini %d >= end %d\n",
      n, i, heap[i].ini, heap[i].end);
  }
}

void heap_swap(HEAP *heap, int i, int j) {

  HEAP temp = heap[i];

  heap[i] = heap[j];

  if (heap[j].prev >= 0) {
    if (heap[j].prev == i) {
      heap[i].prev = j;
    }
    else {
      heap[heap[j].prev].next = i;
    }
  }

  if (heap[j].next >= 0) {
    if (heap[j].next == i) {
      heap[i].next = j;
    }
    else {
      heap[heap[j].next].prev = i;
    }
  }

  heap[j] = temp;
  if (temp.prev >= 0) {
    if (temp.prev == j) {
      heap[j].prev = i;
    }
    else {
      heap[temp.prev].next = j;
    }
  }
  if (temp.next >= 0) {
    if (temp.next == j) {
      heap[j].next = i;
    }
    else {
      heap[temp.next].prev = j;
    }
  }

}

void heap_restore_up(HEAP *heap, int i) {
  int p;
  while (i > 0) {
    p = (i - 2 + (i % 2)) / 2;
    if (heap[i].val < heap[p].val) {
      heap_swap(heap, i, p);
      i = p;
    }
    else {
      break;
    }
  }
}

int heap_select_son (HEAP *heap, int i, int n) {

  int s1 = (2*i)+1;
  int s2 = s1+1;

  if (s2 > n) {
    return -1;
  }
  else if ( (s2 == n) || ( (heap[s1].val) <= (heap[s2].val) ) ) {
    return s1;
  }
  else {
    return s2;
  }
}

void heap_restore_down(HEAP *heap, int i, int n) {
  int s;
  while (true) {
    s = heap_select_son(heap, i, n);
    if ( (s == -1) || (heap[s].val >= heap[i].val) ) {
      break;
    }
    else {
      heap_swap(heap, i, s);
      i = s;
    }
  }
}

void heap_restore(HEAP *heap, int i, int n) {
  int p = (i - 2 + (i % 2)) / 2;
  if ( (p >= 0) && (heap[i].val < heap[p].val) ) {
    heap_restore_up(heap, i);
  }
  else {
    heap_restore_down(heap, i, n);
  }
}

void heap_initialize(HEAP *heap, int n) {
 int i, p, count = 0;
 bool changed = true;
 for (i = n-1; i>= 0; i--) {
   heap_restore_down(heap, i, n);
 }
}

HEAP heap_remove_root(HEAP *heap, int n) {

  HEAP ans = heap[0];
  HEAP new_root = heap[n-1];

  if (ans.next == (n-1) ) {
    ans.next = 0;
  }
  if (ans.prev == (n-1) ) {
    ans.prev = 0;
  }
  if (new_root.prev == 0) {
    new_root.prev = ans.prev;
  }
  if (new_root.next == 0) {
    new_root.next = ans.next;
  }

  heap[0] = new_root;

  if (ans.prev >= 0) {
    heap[ans.prev].next = ans.next;
  }

  if (ans.next >= 0) {
    heap[ans.next].prev = ans.prev;
    heap[ans.next].ini = ans.ini;
    heap[ans.next].ini_vec = ans.ini_vec;
  }

  if (new_root.prev >= 0) {
    heap[new_root.prev].next = 0;
  }
  if (new_root.next >=0) {
    heap[new_root.next].prev = 0;
  }

  heap_restore_down(heap, 0, n-1);

  return ans;
}
