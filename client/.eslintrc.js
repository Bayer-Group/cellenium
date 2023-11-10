module.exports = {
  extends: [
    'airbnb',
    'airbnb-typescript',
    'airbnb/hooks',
    'eslint:recommended',
    'plugin:import/recommended',
    'plugin:react/recommended',
    'plugin:@typescript-eslint/recommended',
    'plugin:prettier/recommended',
  ],
  parserOptions: {
    project: './tsconfig.eslint.json',
  },
  rules: {
    // Disables jsx-a11y https://github.com/import-js/eslint-plugin-import/blob/v2.25.4/docs/rules/no-webpack-loader-syntax.md
    // eslint-disable-next-line global-require
    ...Object.keys(require('eslint-plugin-jsx-a11y').rules).reduce((acc, rule) => {
      acc[`jsx-a11y/${rule}`] = 'off';
      return acc;
    }, {}),
    'react/jsx-uses-react': 'off',
    'react/react-in-jsx-scope': 'off',
    'class-methods-use-this': 'off',
    'linebreak-style': 'off',
    'no-console': 'off',
    'no-continue': 'off',
    'no-multi-assign': 'warn',
    'no-nested-ternary': 'off',
    'no-return-assign': 'warn',
    'no-restricted-exports': 'off',
    'no-restricted-syntax': 'off',
    'no-plusplus': 'off',
    'no-prototype-builtins': 'warn',
    'no-minusminus': 'off',
    'no-underscore-dangle': 'off',
    '@typescript-eslint/no-unused-expressions': [
      'error',
      {
        allowShortCircuit: true,
        allowTernary: true,
        allowTaggedTemplates: true,
      },
    ],
    'max-classes-per-file': 'off',
    'no-param-reassign': ['warn', { props: true, ignorePropertyModificationsFor: ['state'] }], // Exclude state as required by redux-toolkit: https://redux-toolkit.js.org/usage/immer-reducers#linting-state-mutations
    'cypress/unsafe-to-chain-command': 'off',
    'import/no-extraneous-dependencies': 'off',
    'import/no-webpack-loader-syntax': 'off', // Disable to allow webpack file-loaders syntax
    'import/no-unresolved': 'off', // Disable to allow webpack file-loaders syntax
    'import/prefer-default-export': 'off',
    'import/order': 'error',
    'prefer-destructuring': ['warn', { object: true, array: false }],
    'prefer-promise-reject-errors': 'warn',
    'prefer-spread': 'warn',
    '@typescript-eslint/ban-ts-comment': 'warn',
    'react/destructuring-assignment': 'off',
    'react/jsx-curly-brace-presence': 'warn',
    'react/jsx-props-no-spreading': 'off',
    'react/no-unused-class-component-methods': 'warn',
    'react/prop-types': 'off',
    'react/require-default-props': 'off',
    'react/static-property-placement': [
      'warn',
      'property assignment',
      {
        childContextTypes: 'static getter',
        contextTypes: 'static public field',
        contextType: 'static public field',
        displayName: 'static public field',
      },
    ],
  },
  overrides: [
    {
      files: ['{src|tests}/**/*.{test|spec}.ts'],
      extends: ['plugin:jest/recommended'],
      plugins: ['jest'],
    },
  ],
};
