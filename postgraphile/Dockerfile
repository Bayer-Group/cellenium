FROM node:16
WORKDIR /postgraphile
COPY package.json /postgraphile
COPY yarn.lock /postgraphile
RUN yarn
COPY postgraphile_server.js /postgraphile
EXPOSE 6000
ENTRYPOINT ["node", "postgraphile_server.js" , "-n", "0.0.0.0"]
