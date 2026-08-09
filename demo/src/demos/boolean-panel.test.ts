// Tests for the Boolean Gallery sidebar's pure presentation helpers
// (boolean-panel.ts). The DOM assembly itself needs a browser; what is worth
// pinning here is the closed Offset disclosure's summary, which is the only
// place those values are visible once the sliders are folded away.

import { expect, test } from 'bun:test';
import { offsetSummary } from './boolean-panel.ts';

test('the closed offset row still prints its values', () => {
  expect(offsetSummary([0.1, 0, 0])).toEqual({ text: '0.1, 0, 0', muted: false });
  expect(offsetSummary([0.7, -0.2, 0.4])).toEqual({ text: '0.7, -0.2, 0.4', muted: false });
});

test('all-zero offsets read as none', () => {
  expect(offsetSummary([0, 0, 0])).toEqual({ text: 'none', muted: true });
  // -0 is still no offset, but it must not print as "-0".
  expect(offsetSummary([-0, 0, 0]).text).toBe('none');
});
